
package physics3d;

/**
 * Distance joint: keeps two anchor points at a fixed distance.
 *
 * Hot-path optimized: no Vec3/Mat3 allocations in warmStart/solveVelocity.
 */
public final class DistanceJoint implements Joint {
    private final RigidBody A;
    private final RigidBody B;
    private final Vec3 localAnchorA;
    private final Vec3 localAnchorB;

    public float length;
    public float stiffness = 1.0f; // 0..1 (softness; 1 is rigid)

    // Cached (warm start)
    private float impulse = 0f;

    // Constraint direction (from A -> B) in world space
    private final Vec3 n = new Vec3(1, 0, 0);

    // Effective mass along n
    private float effectiveMass = 0f;

    // ---- scratch (no allocations in solver) ----
    private final Vec3 pA = new Vec3();
    private final Vec3 pB = new Vec3();
    private final Vec3 d = new Vec3();

    private final Vec3 rA = new Vec3();
    private final Vec3 rB = new Vec3();

    private final Vec3 rAxn = new Vec3();
    private final Vec3 rBxn = new Vec3();

    private final Vec3 invI_rAxn = new Vec3();
    private final Vec3 invI_rBxn = new Vec3();
    private final Vec3 tmpCross = new Vec3();

    private final Vec3 wXrA = new Vec3();
    private final Vec3 wXrB = new Vec3();
    private final Vec3 vA = new Vec3();
    private final Vec3 vB = new Vec3();
    private final Vec3 rv = new Vec3();

    private final Vec3 P = new Vec3();
    private final Vec3 negP = new Vec3();

    // inv inertia helper scratch
    private final Vec3 body = new Vec3();
    private final Vec3 scaled = new Vec3();

    public DistanceJoint(RigidBody a, RigidBody b, Vec3 worldAnchorA, Vec3 worldAnchorB) {
        this.A = a;
        this.B = b;

        // constructor allocations ok
        this.localAnchorA = a.localPointFromWorld(worldAnchorA);
        this.localAnchorB = b.localPointFromWorld(worldAnchorB);

        // initial length
        this.length = worldAnchorB.sub(worldAnchorA).len();
    }

    @Override public RigidBody bodyA() { return A; }
    @Override public RigidBody bodyB() { return B; }

    /**
     * Computes current n and effective mass (allocation-free).
     */
    private void precompute() {
        // world anchors
        A.worldPointFromLocalInto(pA, localAnchorA);
        B.worldPointFromLocalInto(pB, localAnchorB);

        // d = pB - pA
        d.set(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
        float dist2 = d.len2();
        float dist = (float) Math.sqrt(dist2);

        if (dist < 1e-6f) {
            n.set(1, 0, 0);
        } else {
            float inv = 1f / dist;
            n.set(d.x * inv, d.y * inv, d.z * inv);
        }

        // rA = pA - A.pos, rB = pB - B.pos
        rA.set(pA.x - A.position.x, pA.y - A.position.y, pA.z - A.position.z);
        rB.set(pB.x - B.position.x, pB.y - B.position.y, pB.z - B.position.z);

        // rAxn = rA x n, rBxn = rB x n
        rAxn.setCross(rA, n);
        rBxn.setCross(rB, n);

        // k = invMassA + invMassB + n · ( (invI*(rAxn))×rA + (invI*(rBxn))×rB )
        float k = A.invMass + B.invMass;

        invInertiaWorldMulInto(invI_rAxn, A, rAxn);
        tmpCross.setCross(invI_rAxn, rA);
        k += n.dot(tmpCross);

        invInertiaWorldMulInto(invI_rBxn, B, rBxn);
        tmpCross.setCross(invI_rBxn, rB);
        k += n.dot(tmpCross);

        effectiveMass = 1f / Math.max(1e-8f, k);
    }

    @Override
    public void warmStart() {
        if (A.isStatic() && B.isStatic()) return;

        precompute();

        // Apply cached impulse along n at anchors
        P.set(n.x * impulse, n.y * impulse, n.z * impulse);
        negP.set(-P.x, -P.y, -P.z);

        A.applyImpulse(negP, pA);
        B.applyImpulse(P, pB);
    }

    @Override
    public void solveVelocity(float dt) {
        if (A.isStatic() && B.isStatic()) return;

        precompute();

        // Relative velocity at anchors:
        // vA = vA + wA×rA, vB = vB + wB×rB
        wXrA.setCross(A.angularVelocity, rA);
        vA.set(A.linearVelocity.x + wXrA.x, A.linearVelocity.y + wXrA.y, A.linearVelocity.z + wXrA.z);

        wXrB.setCross(B.angularVelocity, rB);
        vB.set(B.linearVelocity.x + wXrB.x, B.linearVelocity.y + wXrB.y, B.linearVelocity.z + wXrB.z);

        rv.set(vB.x - vA.x, rv.y = vB.y - vA.y, rv.z = vB.z - vA.z);

        float Cdot = rv.dot(n);

        // positional error as spring bias (soft)
        float dist = (float) Math.sqrt(d.len2());
        float C = dist - length;
        float bias = -stiffness * (C / Math.max(1e-6f, dt));

        float lambda = -(Cdot + bias) * effectiveMass;
        impulse += lambda;

        // Apply impulse
        P.set(n.x * lambda, n.y * lambda, n.z * lambda);
        negP.set(-P.x, -P.y, -P.z);

        A.applyImpulse(negP, pA);
        B.applyImpulse(P, pB);
    }

    @Override
    public void solvePosition() {
        // optional: mild positional correction (kept empty for now)
    }

    // --------------------------------------------------------------------
    // invI_world * v (allocation-free): R( invI_body * (R^T v) ), invI_body diagonal
    // --------------------------------------------------------------------
    private void invInertiaWorldMulInto(Vec3 out, RigidBody rb, Vec3 vWorld) {
        rotateConjugateInto(body, rb.orientation, vWorld);

        // invInertiaBody is diagonal in this engine
        scaled.set(body.x * rb.invInertiaBody.m00, body.y * rb.invInertiaBody.m11, body.z * rb.invInertiaBody.m22);

        rotateInto(out, rb.orientation, scaled);
    }

    private static void rotateInto(Vec3 out, Quat q, Vec3 v) {
        final float qw = q.w, qx = q.x, qy = q.y, qz = q.z;
        final float vx = v.x, vy = v.y, vz = v.z;

        final float tx = 2f * (qy * vz - qz * vy);
        final float ty = 2f * (qz * vx - qx * vz);
        final float tz = 2f * (qx * vy - qy * vx);

        out.set(
                vx + qw * tx + (qy * tz - qz * ty),
                vy + qw * ty + (qz * tx - qx * tz),
                vz + qw * tz + (qx * ty - qy * tx)
        );
    }

    private static void rotateConjugateInto(Vec3 out, Quat q, Vec3 v) {
        final float qw = q.w, qx = q.x, qy = q.y, qz = q.z;
        final float vx = v.x, vy = v.y, vz = v.z;

        final float tx = 2f * (qy * vz - qz * vy);
        final float ty = 2f * (qz * vx - qx * vz);
        final float tz = 2f * (qx * vy - qy * vx);

        out.set(
                vx - qw * tx + (qy * tz - qz * ty),
                vy - qw * ty + (qz * tx - qx * tz),
                vz - qw * tz + (qx * ty - qy * tx)
        );
    }
}
