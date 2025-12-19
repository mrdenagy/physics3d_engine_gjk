
package physics3d;

/**
 * Slider (prismatic) joint:
 * - Perpendicular translation constrained with a full 3x3 effective mass solve,
 *   then projected to remove the allowed axis component (so motion along axis remains).
 * - Optional translation limits along axis.
 * - Optional motor along axis.
 *
 * Hot-path optimized: no Vec3/Mat3 allocations during warmStart/solveVelocity.
 */
public final class SliderJoint implements Joint {
    private final RigidBody A;
    private final RigidBody B;
    private final Vec3 localAnchorA;
    private final Vec3 localAnchorB;
    private final Vec3 localAxisA;

    // cached impulses
    private final Vec3 impulsePerp = new Vec3(0, 0, 0);
    private float motorImpulse = 0f;
    private float lowerImpulse = 0f;
    private float upperImpulse = 0f;

    public float softness = 0.0f;

    // Limits and motor
    public boolean enableLimit = false;
    public float lowerTranslation = -1.0f;
    public float upperTranslation = 1.0f;

    public boolean enableMotor = false;
    public float motorSpeed = 0f;     // m/s along axis
    public float maxMotorForce = 10f; // N

    // ---- scratch (no allocations in solver) ----
    private final Vec3 ax = new Vec3();
    private final Vec3 pA = new Vec3();
    private final Vec3 pB = new Vec3();
    private final Vec3 rA = new Vec3();
    private final Vec3 rB = new Vec3();

    private final Vec3 wXrA = new Vec3();
    private final Vec3 wXrB = new Vec3();
    private final Vec3 vA = new Vec3();
    private final Vec3 vB = new Vec3();
    private final Vec3 Cdot = new Vec3();

    private final Vec3 lambda = new Vec3();
    private final Vec3 negLambda = new Vec3();

    private final Vec3 imp = new Vec3();
    private final Vec3 negImp = new Vec3();

    private final Vec3 d = new Vec3();       // pB - pA
    private final Vec3 tmp = new Vec3();     // general scratch
    private final Vec3 body = new Vec3();    // inv inertia helper
    private final Vec3 scaled = new Vec3();  // inv inertia helper

    public SliderJoint(RigidBody a, RigidBody b, Vec3 worldAnchor, Vec3 worldAxis) {
        this.A = a;
        this.B = b;

        // one-time allocations ok in constructor
        this.localAnchorA = a.localPointFromWorld(worldAnchor);
        this.localAnchorB = b.localPointFromWorld(worldAnchor);

        Vec3 axisW = worldAxis.normalized();
        // store axis in A's local space (direction only)
        Vec3 p1 = a.localPointFromWorld(a.position.add(axisW));
        Vec3 p0 = a.localPointFromWorld(a.position);
        this.localAxisA = p1.sub(p0).normalized();
    }

    @Override public RigidBody bodyA() { return A; }
    @Override public RigidBody bodyB() { return B; }

    // ---- allocation-free helpers ----

    private void axisWorldInto(Vec3 out) {
        // out = normalize( A.orientation * localAxisA )
        rotateInto(out, A.orientation, localAxisA);
        normalizeSelf(out);
    }

    private void anchorWorldAInto(Vec3 out) { A.worldPointFromLocalInto(out, localAnchorA); }
    private void anchorWorldBInto(Vec3 out) { B.worldPointFromLocalInto(out, localAnchorB); }

    @Override
    public void warmStart() {
        if (A.isStatic() && B.isStatic()) return;

        anchorWorldAInto(pA);
        anchorWorldBInto(pB);

        // apply cached perpendicular impulse at anchors
        negImp.set(-impulsePerp.x, -impulsePerp.y, -impulsePerp.z);
        A.applyImpulse(negImp, pA);
        B.applyImpulse(impulsePerp, pB);

        axisWorldInto(ax);
        float sum = motorImpulse + lowerImpulse + upperImpulse;
        if (sum != 0f) {
            imp.set(ax.x * sum, ax.y * sum, ax.z * sum);
            negImp.set(-imp.x, -imp.y, -imp.z);
            A.applyImpulse(negImp, pA);
            B.applyImpulse(imp, pB);
        }
    }

    @Override
    public void solveVelocity(float dt) {
        if (A.isStatic() && B.isStatic()) return;

        axisWorldInto(ax);
        anchorWorldAInto(pA);
        anchorWorldBInto(pB);

        rA.set(pA.x - A.position.x, pA.y - A.position.y, pA.z - A.position.z);
        rB.set(pB.x - B.position.x, pB.y - B.position.y, pB.z - B.position.z);

        // Cdot = vB + wB×rB - (vA + wA×rA)
        wXrA.setCross(A.angularVelocity, rA);
        vA.set(A.linearVelocity.x + wXrA.x, A.linearVelocity.y + wXrA.y, A.linearVelocity.z + wXrA.z);

        wXrB.setCross(B.angularVelocity, rB);
        vB.set(B.linearVelocity.x + wXrB.x, B.linearVelocity.y + wXrB.y, B.linearVelocity.z + wXrB.z);

        Cdot.set(vB.x - vA.x, vB.y - vA.y, vB.z - vA.z);

        // Solve 3x3 point-to-point effective mass for lambda = -K^-1 * Cdot
        solvePointToPointLambda(lambda, A, rA, B, rB, Cdot);

        // Project out allowed axis component so perpendicular is constrained
        float along = lambda.dot(ax);
        lambda.set(lambda.x - ax.x * along, lambda.y - ax.y * along, lambda.z - ax.z * along);

        lambda.mulLocal(1f - softness);

        impulsePerp.addLocal(lambda);

        negLambda.set(-lambda.x, -lambda.y, -lambda.z);
        A.applyImpulse(negLambda, pA);
        B.applyImpulse(lambda, pB);

        // ---- Motor along axis ----
        if (enableMotor) {
            // relative speed along axis
            wXrA.setCross(A.angularVelocity, rA);
            vA.set(A.linearVelocity.x + wXrA.x, A.linearVelocity.y + wXrA.y, A.linearVelocity.z + wXrA.z);

            wXrB.setCross(B.angularVelocity, rB);
            vB.set(B.linearVelocity.x + wXrB.x, B.linearVelocity.y + wXrB.y, B.linearVelocity.z + wXrB.z);

            Cdot.set(vB.x - vA.x, vB.y - vA.y, vB.z - vA.z);
            float speed = Cdot.dot(ax);

            float CdotMotor = speed - motorSpeed;
            float kAxis = axisEffectiveMass(A, rA, ax) + axisEffectiveMass(B, rB, ax);
            float mEff = 1f / Math.max(1e-8f, kAxis);

            float lam = -CdotMotor * mEff;
            float old = motorImpulse;
            float maxImp = maxMotorForce * dt;
            motorImpulse = MathUtil.clamp(old + lam, -maxImp, maxImp);
            float dImp = motorImpulse - old;

            imp.set(ax.x * dImp, ax.y * dImp, ax.z * dImp);
            negImp.set(-imp.x, -imp.y, -imp.z);
            A.applyImpulse(negImp, pA);
            B.applyImpulse(imp, pB);
        }

        // ---- Limits along axis ----
        if (enableLimit) {
            float trans = currentTranslation(ax);

            // lower
            if (trans < lowerTranslation) {
                float C = trans - lowerTranslation;
                float bias = -0.2f * (C / Math.max(1e-6f, dt));

                // recompute rel speed
                wXrA.setCross(A.angularVelocity, rA);
                vA.set(A.linearVelocity.x + wXrA.x, A.linearVelocity.y + wXrA.y, A.linearVelocity.z + wXrA.z);

                wXrB.setCross(B.angularVelocity, rB);
                vB.set(B.linearVelocity.x + wXrB.x, B.linearVelocity.y + wXrB.y, B.linearVelocity.z + wXrB.z);

                Cdot.set(vB.x - vA.x, vB.y - vA.y, vB.z - vA.z);
                float rel = Cdot.dot(ax);

                float kAxis = axisEffectiveMass(A, rA, ax) + axisEffectiveMass(B, rB, ax);
                float mEff = 1f / Math.max(1e-8f, kAxis);

                float lam = -(rel + bias) * mEff;

                float old = lowerImpulse;
                lowerImpulse = Math.max(old + lam, 0f);
                float dL = lowerImpulse - old;

                imp.set(ax.x * dL, ax.y * dL, ax.z * dL);
                negImp.set(-imp.x, -imp.y, -imp.z);
                A.applyImpulse(negImp, pA);
                B.applyImpulse(imp, pB);
            } else {
                lowerImpulse *= 0.5f;
            }

            // upper
            if (trans > upperTranslation) {
                float C = trans - upperTranslation;
                float bias = -0.2f * (C / Math.max(1e-6f, dt));

                wXrA.setCross(A.angularVelocity, rA);
                vA.set(A.linearVelocity.x + wXrA.x, A.linearVelocity.y + wXrA.y, A.linearVelocity.z + wXrA.z);

                wXrB.setCross(B.angularVelocity, rB);
                vB.set(B.linearVelocity.x + wXrB.x, B.linearVelocity.y + wXrB.y, B.linearVelocity.z + wXrB.z);

                Cdot.set(vB.x - vA.x, vB.y - vA.y, vB.z - vA.z);
                float rel = Cdot.dot(ax);

                float kAxis = axisEffectiveMass(A, rA, ax) + axisEffectiveMass(B, rB, ax);
                float mEff = 1f / Math.max(1e-8f, kAxis);

                float lam = -(rel + bias) * mEff;

                float old = upperImpulse;
                upperImpulse = Math.min(old + lam, 0f);
                float dU = upperImpulse - old;

                imp.set(ax.x * dU, ax.y * dU, ax.z * dU);
                negImp.set(-imp.x, -imp.y, -imp.z);
                A.applyImpulse(negImp, pA);
                B.applyImpulse(imp, pB);
            } else {
                upperImpulse *= 0.5f;
            }
        }
    }

    private float currentTranslation(Vec3 axWorld) {
        anchorWorldAInto(pA);
        anchorWorldBInto(pB);
        d.set(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
        return d.dot(axWorld);
    }

    // -------- Effective mass helpers (allocation-free) --------

    private void solvePointToPointLambda(Vec3 outLambda, RigidBody A, Vec3 rA, RigidBody B, Vec3 rB, Vec3 Cdot) {
        // K = (imA+imB)I - [rA] invI_A [rA]^T - [rB] invI_B [rB]^T
        float im = A.invMass + B.invMass;

        // Build angular term for a body: T = [r] invI [r]^T, where invI_world applied via helper
        // We compute columns by applying T to basis vectors e1,e2,e3 via:
        // T * x = (invI*(r×x)) × r
        float K00 = im, K01 = 0,  K02 = 0;
        float K10 = 0,  K11 = im, K12 = 0;
        float K20 = 0,  K21 = 0,  K22 = im;

        // Add body A contribution
        addAngularK(A, rA, 1, 0, 0); K00 -= tmp.x; K10 -= tmp.y; K20 -= tmp.z;
        addAngularK(A, rA, 0, 1, 0); K01 -= tmp.x; K11 -= tmp.y; K21 -= tmp.z;
        addAngularK(A, rA, 0, 0, 1); K02 -= tmp.x; K12 -= tmp.y; K22 -= tmp.z;

        // Add body B contribution
        addAngularK(B, rB, 1, 0, 0); K00 -= tmp.x; K10 -= tmp.y; K20 -= tmp.z;
        addAngularK(B, rB, 0, 1, 0); K01 -= tmp.x; K11 -= tmp.y; K21 -= tmp.z;
        addAngularK(B, rB, 0, 0, 1); K02 -= tmp.x; K12 -= tmp.y; K22 -= tmp.z;

        // Solve outLambda = -inv(K) * Cdot
        float bx = -Cdot.x, by = -Cdot.y, bz = -Cdot.z;
        solve3x3(outLambda, K00, K01, K02, K10, K11, K12, K20, K21, K22, bx, by, bz);
    }

    // tmp = (invI_world*(r×basis))×r
    private void addAngularK(RigidBody bodyRB, Vec3 r, float bx, float by, float bz) {
        // r×basis
        tmp.setCross(r, tmp.set(bx, by, bz));
        // invI_world * (r×basis)
        invInertiaWorldMulInto(body, bodyRB, tmp);
        // (invI*(...)) × r
        tmp.setCross(body, r);
    }

    private float axisEffectiveMass(RigidBody bodyRB, Vec3 r, Vec3 axis) {
        float im = bodyRB.invMass;
        // ang term = axis · ( (invI*(r×axis)) × r )
        tmp.setCross(r, axis);
        invInertiaWorldMulInto(body, bodyRB, tmp);
        tmp.setCross(body, r);
        float ang = axis.dot(tmp);
        return im + ang;
    }

    // invI_world * v = R( invI_body * (R^T v) ), invI_body diagonal
    private void invInertiaWorldMulInto(Vec3 out, RigidBody rb, Vec3 vWorld) {
        rotateConjugateInto(body, rb.orientation, vWorld);
        scaled.set(body.x * rb.invInertiaBody.m00, body.y * rb.invInertiaBody.m11, body.z * rb.invInertiaBody.m22);
        rotateInto(out, rb.orientation, scaled);
    }

    // Solve 3x3 linear system (general) using explicit inverse (Cramer's / adjugate)
    private static void solve3x3(Vec3 out,
                                 float a00, float a01, float a02,
                                 float a10, float a11, float a12,
                                 float a20, float a21, float a22,
                                 float bx, float by, float bz) {

        float c00 = a11 * a22 - a12 * a21;
        float c01 = a02 * a21 - a01 * a22;
        float c02 = a01 * a12 - a02 * a11;

        float c10 = a12 * a20 - a10 * a22;
        float c11 = a00 * a22 - a02 * a20;
        float c12 = a02 * a10 - a00 * a12;

        float c20 = a10 * a21 - a11 * a20;
        float c21 = a01 * a20 - a00 * a21;
        float c22 = a00 * a11 - a01 * a10;

        float det = a00 * c00 + a01 * c10 + a02 * c20;
        float invDet = (Math.abs(det) < 1e-10f) ? 0f : (1f / det);

        out.set(
                (c00 * bx + c01 * by + c02 * bz) * invDet,
                (c10 * bx + c11 * by + c12 * bz) * invDet,
                (c20 * bx + c21 * by + c22 * bz) * invDet
        );
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

    private static void normalizeSelf(Vec3 v) {
        float l2 = v.x * v.x + v.y * v.y + v.z * v.z;
        if (l2 < 1e-12f) return;
        float inv = 1f / (float) Math.sqrt(l2);
        v.set(v.x * inv, v.y * inv, v.z * inv);
    }

    @Override
    public void solvePosition() {
        // optional
    }
}
