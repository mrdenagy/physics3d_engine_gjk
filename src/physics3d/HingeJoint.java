
package physics3d;

/**
 * Hinge (revolute) joint:
 * - Coincident anchors (3 DOF) solved with full 3x3 effective mass (explicit float inverse).
 * - Axis alignment (2 DOF) solved by angular impulses about two orthogonal axes.
 * - Optional angular motor and limits about hinge axis.
 *
 * Hot-path optimized: no Vec3/Mat3 allocations during warmStart/solveVelocity.
 */
public final class HingeJoint implements Joint {
    private final RigidBody A;
    private final RigidBody B;

    private final Vec3 localAnchorA;
    private final Vec3 localAnchorB;

    private final Vec3 localAxisA;
    private final Vec3 localAxisB;

    // Cached impulses
    private final Vec3 impulseP2P = new Vec3(0, 0, 0);
    private final Vec3 impulseAlign = new Vec3(0, 0, 0); // world-space angular impulse (perp to hinge axis)

    private float motorImpulse = 0f;
    private float limitImpulse = 0f;

    // Motor/limit settings
    public boolean enableMotor = false;
    public float motorSpeed = 0f; // rad/s about hinge axis
    public float maxMotorTorque = 5f; // N*m

    public boolean enableLimit = false;
    public float lowerAngle = (float)(-Math.PI/4);
    public float upperAngle = (float)( Math.PI/4);

    public float softness = 0.0f;

    // ---- scratch ----
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

    private final Vec3 hA = new Vec3();
    private final Vec3 hB = new Vec3();
    private final Vec3 h = new Vec3();

    private final Vec3 cross = new Vec3();
    private final Vec3 u = new Vec3();
    private final Vec3 v = new Vec3();

    private final Vec3 wRel = new Vec3();

    private final Vec3 angImp = new Vec3();
    private final Vec3 negAngImp = new Vec3();

    // inv inertia helper scratch
    private final Vec3 body = new Vec3();
    private final Vec3 scaled = new Vec3();
    private final Vec3 tmp = new Vec3();

    public HingeJoint(RigidBody a, RigidBody b, Vec3 worldAnchor, Vec3 worldAxis) {
        this.A = a;
        this.B = b;

        this.localAnchorA = a.localPointFromWorld(worldAnchor);
        this.localAnchorB = b.localPointFromWorld(worldAnchor);

        Vec3 axisW = worldAxis.normalized();

        // represent local axis by transforming two points (constructor allocations ok)
        this.localAxisA = a.localPointFromWorld(a.position.add(axisW)).sub(a.localPointFromWorld(a.position)).normalized();
        this.localAxisB = b.localPointFromWorld(b.position.add(axisW)).sub(b.localPointFromWorld(b.position)).normalized();
    }

    @Override public RigidBody bodyA() { return A; }
    @Override public RigidBody bodyB() { return B; }

    private void axisWorldAInto(Vec3 out) { rotateInto(out, A.orientation, localAxisA); normalizeSelf(out); }
    private void axisWorldBInto(Vec3 out) { rotateInto(out, B.orientation, localAxisB); normalizeSelf(out); }

    private void anchorWorldAInto(Vec3 out) { A.worldPointFromLocalInto(out, localAnchorA); }
    private void anchorWorldBInto(Vec3 out) { B.worldPointFromLocalInto(out, localAnchorB); }

    @Override
    public void warmStart() {
        if (A.isStatic() && B.isStatic()) return;

        anchorWorldAInto(pA);
        anchorWorldBInto(pB);

        // point-to-point cached impulse
        negLambda.set(-impulseP2P.x, -impulseP2P.y, -impulseP2P.z);
        A.applyImpulse(negLambda, pA);
        B.applyImpulse(impulseP2P, pB);

        // axis alignment cached angular impulse (apply directly to angular velocity)
        if (!A.isStatic()) {
            invInertiaWorldMulInto(tmp, A, impulseAlign);
            tmp.negLocal();
            A.angularVelocity.addLocal(tmp);
        }
        if (!B.isStatic()) {
            invInertiaWorldMulInto(tmp, B, impulseAlign);
            B.angularVelocity.addLocal(tmp);
        }

        // motor + limit cached about hinge axis
        axisWorldAInto(h);
        float sum = motorImpulse + limitImpulse;
        if (sum != 0f) {
            angImp.set(h.x * sum, h.y * sum, h.z * sum);

            if (!A.isStatic()) {
                invInertiaWorldMulInto(tmp, A, angImp);
                tmp.negLocal();
                A.angularVelocity.addLocal(tmp);
            }
            if (!B.isStatic()) {
                invInertiaWorldMulInto(tmp, B, angImp);
                B.angularVelocity.addLocal(tmp);
            }
        }
    }

    @Override
    public void solveVelocity(float dt) {
        if (A.isStatic() && B.isStatic()) return;

        anchorWorldAInto(pA);
        anchorWorldBInto(pB);

        rA.set(pA.x - A.position.x, pA.y - A.position.y, pA.z - A.position.z);
        rB.set(pB.x - B.position.x, pB.y - B.position.y, pB.z - B.position.z);

        // --- Point-to-point (3DOF) ---
        wXrA.setCross(A.angularVelocity, rA);
        vA.set(A.linearVelocity.x + wXrA.x, A.linearVelocity.y + wXrA.y, A.linearVelocity.z + wXrA.z);

        wXrB.setCross(B.angularVelocity, rB);
        vB.set(B.linearVelocity.x + wXrB.x, B.linearVelocity.y + wXrB.y, B.linearVelocity.z + wXrB.z);

        Cdot.set(vB.x - vA.x, vB.y - vA.y, vB.z - vA.z);

        solvePointToPointLambda(lambda, A, rA, B, rB, Cdot);
        lambda.mulLocal(1f - softness);

        impulseP2P.addLocal(lambda);

        negLambda.set(-lambda.x, -lambda.y, -lambda.z);
        A.applyImpulse(negLambda, pA);
        B.applyImpulse(lambda, pB);

        // --- Axis alignment (2DOF) ---
        axisWorldAInto(hA);
        axisWorldBInto(hB);

        cross.setCross(hA, hB);
        float crossLen2 = cross.len2();
        if (crossLen2 > 1e-12f) {
            float crossLen = (float)Math.sqrt(crossLen2);
            u.set(cross.x / crossLen, cross.y / crossLen, cross.z / crossLen); // perp to both
            v.setCross(hA, u);
            normalizeSelf(v);

            wRel.set(B.angularVelocity.x - A.angularVelocity.x,
                     B.angularVelocity.y - A.angularVelocity.y,
                     B.angularVelocity.z - A.angularVelocity.z);

            solveAngular1D(u, wRel, dt);
            solveAngular1D(v, wRel, dt);
        }

        // --- Motor and limits about hinge axis ---
        h.set(hA);

        wRel.set(B.angularVelocity.x - A.angularVelocity.x,
                 B.angularVelocity.y - A.angularVelocity.y,
                 B.angularVelocity.z - A.angularVelocity.z);
        float wRelH = wRel.dot(h);

        // motor
        if (enableMotor) {
            float CdotMotor = wRelH - motorSpeed;
            float k = axisAngularMass(A, h) + axisAngularMass(B, h);
            float mEff = 1f / Math.max(1e-8f, k);

            float lam = -CdotMotor * mEff;

            float old = motorImpulse;
            float maxImp = maxMotorTorque * dt;
            motorImpulse = MathUtil.clamp(old + lam, -maxImp, maxImp);
            float dImp = motorImpulse - old;

            angImp.set(h.x * dImp, h.y * dImp, h.z * dImp);

            if (!A.isStatic()) {
                invInertiaWorldMulInto(tmp, A, angImp);
                tmp.negLocal();
                A.angularVelocity.addLocal(tmp);
            }
            if (!B.isStatic()) {
                invInertiaWorldMulInto(tmp, B, angImp);
                B.angularVelocity.addLocal(tmp);
            }
        }

        // limit
        if (enableLimit) {
            float angle = currentAngle(h);
            float C = 0f;
            int limitState = 0; // -1 lower, +1 upper

            if (angle < lowerAngle) { C = angle - lowerAngle; limitState = -1; }
            else if (angle > upperAngle) { C = angle - upperAngle; limitState = 1; }

            if (limitState != 0) {
                float bias = -0.2f * (C / Math.max(1e-6f, dt));
                float k = axisAngularMass(A, h) + axisAngularMass(B, h);
                float mEff = 1f / Math.max(1e-8f, k);

                float lam = -(wRelH + bias) * mEff;

                float old = limitImpulse;
                if (limitState < 0) limitImpulse = Math.max(old + lam, 0f);
                else limitImpulse = Math.min(old + lam, 0f);

                float dImp = limitImpulse - old;

                angImp.set(h.x * dImp, h.y * dImp, h.z * dImp);

                if (!A.isStatic()) {
                    invInertiaWorldMulInto(tmp, A, angImp);
                    tmp.negLocal();
                    A.angularVelocity.addLocal(tmp);
                }
                if (!B.isStatic()) {
                    invInertiaWorldMulInto(tmp, B, angImp);
                    B.angularVelocity.addLocal(tmp);
                }
            } else {
                limitImpulse *= 0.5f;
            }
        }
    }

    private void solveAngular1D(Vec3 axis, Vec3 wRel, float dt) {
        float Cdot1 = wRel.dot(axis);
        float k = axisAngularMass(A, axis) + axisAngularMass(B, axis);
        float mEff = 1f / Math.max(1e-8f, k);

        float lam = -Cdot1 * mEff;
        lam *= (1f - softness);

        // angular impulse = axis * lam
        angImp.set(axis.x * lam, axis.y * lam, axis.z * lam);
        impulseAlign.addLocal(angImp);

        if (!A.isStatic()) {
            invInertiaWorldMulInto(tmp, A, angImp);
            tmp.negLocal();
            A.angularVelocity.addLocal(tmp);
        }
        if (!B.isStatic()) {
            invInertiaWorldMulInto(tmp, B, angImp);
            B.angularVelocity.addLocal(tmp);
        }
    }

    private float axisAngularMass(RigidBody rb, Vec3 axis) {
        // axis · (invI_world * axis)
        invInertiaWorldMulInto(tmp, rb, axis);
        return axis.dot(tmp);
    }

    private float currentAngle(Vec3 hingeAxisWorld) {
        Quat qRel = Quat.relative(A.orientation, B.orientation);
        return Quat.twistAngleAroundAxis(qRel, hingeAxisWorld);
    }

    // 3x3 anchor effective mass solve (same technique used in SliderJoint)
    private void solvePointToPointLambda(Vec3 outLambda, RigidBody A, Vec3 rA, RigidBody B, Vec3 rB, Vec3 Cdot) {
        float im = A.invMass + B.invMass;

        float K00 = im, K01 = 0,  K02 = 0;
        float K10 = 0,  K11 = im, K12 = 0;
        float K20 = 0,  K21 = 0,  K22 = im;

        addAngularK(A, rA, 1, 0, 0); K00 -= tmp.x; K10 -= tmp.y; K20 -= tmp.z;
        addAngularK(A, rA, 0, 1, 0); K01 -= tmp.x; K11 -= tmp.y; K21 -= tmp.z;
        addAngularK(A, rA, 0, 0, 1); K02 -= tmp.x; K12 -= tmp.y; K22 -= tmp.z;

        addAngularK(B, rB, 1, 0, 0); K00 -= tmp.x; K10 -= tmp.y; K20 -= tmp.z;
        addAngularK(B, rB, 0, 1, 0); K01 -= tmp.x; K11 -= tmp.y; K21 -= tmp.z;
        addAngularK(B, rB, 0, 0, 1); K02 -= tmp.x; K12 -= tmp.y; K22 -= tmp.z;

        float bx = -Cdot.x, by = -Cdot.y, bz = -Cdot.z;
        solve3x3(outLambda, K00, K01, K02, K10, K11, K12, K20, K21, K22, bx, by, bz);
    }

    private void addAngularK(RigidBody rb, Vec3 r, float bx, float by, float bz) {
        // tmp = (invI*(r×basis))×r
        tmp.setCross(r, tmp.set(bx, by, bz));
        invInertiaWorldMulInto(body, rb, tmp);
        tmp.setCross(body, r);
    }

    private void invInertiaWorldMulInto(Vec3 out, RigidBody rb, Vec3 vWorld) {
        rotateConjugateInto(body, rb.orientation, vWorld);
        scaled.set(body.x * rb.invInertiaBody.m00, body.y * rb.invInertiaBody.m11, body.z * rb.invInertiaBody.m22);
        rotateInto(out, rb.orientation, scaled);
    }

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
        // Optional: keep relying on velocity solve.
    }
}
