package physics3d;

/**
 * Continuous collision detection (CCD) using conservative advancement with GJK distance.
 *
 * Rotation-aware conservative advancement and optional angular-motion bound.
 * Allocation-reduced: avoids per-iteration Vec3 allocations.
 */
public final class CCD {
    private CCD(){}

    public static final class TOIResult {
        public boolean hit;
        public float t; // time in [0, dt]
        public Vec3 normal = new Vec3(0,1,0); // from A to B
        public float distance;
    }

    public static TOIResult timeOfImpact(RigidBody A, RigidBody B, float dt, boolean useAngularBound){
        TOIResult out = new TOIResult();
        out.hit = false;
        out.t = dt;

        float vrx = A.linearVelocity.x - B.linearVelocity.x;
        float vry = A.linearVelocity.y - B.linearVelocity.y;
        float vrz = A.linearVelocity.z - B.linearVelocity.z;
        float vRel2 = vrx*vrx + vry*vry + vrz*vrz;

        float wRel2 = A.angularVelocity.len2() + B.angularVelocity.len2();
        if (vRel2 < 1e-12f && wRel2 < 1e-12f) return out;

        Vec3 pA0 = new Vec3(A.position.x, A.position.y, A.position.z);
        Vec3 pB0 = new Vec3(B.position.x, B.position.y, B.position.z);
        Quat qA0 = new Quat(A.orientation.w, A.orientation.x, A.orientation.y, A.orientation.z);
        Quat qB0 = new Quat(B.orientation.w, B.orientation.x, B.orientation.y, B.orientation.z);

        float t = 0f;
        float target = A.shape.margin() + B.shape.margin();
        float tol = 1e-3f;

        float RA = boundingRadius(A);
        float RB = boundingRadius(B);

        Vec3 n = new Vec3(1,0,0);

        for (int iter=0; iter<30; iter++){
            A.position.set(
                    pA0.x + A.linearVelocity.x * t,
                    pA0.y + A.linearVelocity.y * t,
                    pA0.z + A.linearVelocity.z * t
            );
            B.position.set(
                    pB0.x + B.linearVelocity.x * t,
                    pB0.y + B.linearVelocity.y * t,
                    pB0.z + B.linearVelocity.z * t
            );

            A.orientation = quatAt(qA0, A.angularVelocity, t);
            B.orientation = quatAt(qB0, B.angularVelocity, t);

            GJKDistance.Result dres = GJKDistance.distance(A, B);

            if (dres.intersect || dres.distance <= target + tol){
                out.hit = true;
                out.t = t;

                float nx = dres.normal.x, ny = dres.normal.y, nz = dres.normal.z;
                float nl2 = nx*nx + ny*ny + nz*nz;
                if (nl2 > 1e-10f){
                    float inv = 1f/(float)Math.sqrt(nl2);
                    out.normal.set(nx*inv, ny*inv, nz*inv);
                } else {
                    out.normal.set(n);
                }

                out.distance = dres.distance;
                break;
            }

            float dx = dres.normal.x, dy = dres.normal.y, dz = dres.normal.z;
            float dl2 = dx*dx + dy*dy + dz*dz;

            if (dl2 < 1e-10f){
                dx = dres.pointB.x - dres.pointA.x;
                dy = dres.pointB.y - dres.pointA.y;
                dz = dres.pointB.z - dres.pointA.z;
                dl2 = dx*dx + dy*dy + dz*dz;
            }

            if (dl2 < 1e-10f){
                dx = 1f; dy = 0f; dz = 0f;
                dl2 = 1f;
            }

            float invD = 1f/(float)Math.sqrt(dl2);
            n.set(dx*invD, dy*invD, dz*invD);

            float closing = -(vrx*n.x + vry*n.y + vrz*n.z);
            if (useAngularBound){
                float wBound = A.angularVelocity.len()*RA + B.angularVelocity.len()*RB;
                closing += wBound;
            }

            if (closing <= 1e-6f) break;

            float step = (dres.distance - target) / closing;
            step = MathUtil.clamp(step, 0.0f, dt - t);
            if (step < 1e-5f) step = 1e-5f;

            t += step;
            if (t > dt) break;
        }

        A.position.set(pA0);
        B.position.set(pB0);
        A.orientation = qA0;
        B.orientation = qB0;

        return out;
    }

    public static TOIResult timeOfImpact(RigidBody A, RigidBody B, float dt){
        return timeOfImpact(A, B, dt, true);
    }

    private static Quat quatAt(Quat q0, Vec3 omegaWorld, float t){
        if (omegaWorld.len2() < 1e-12f || t == 0f) return new Quat(q0.w, q0.x, q0.y, q0.z);
        Quat dq = Quat.fromAngularVelocity(omegaWorld, t);
        return dq.mul(q0).normalized();
    }

    private static float boundingRadius(RigidBody b){
        if (b.shape.type() == Shape.Type.SPHERE){
            return ((SphereShape)b.shape).r;
        }
        if (b.shape.type() == Shape.Type.BOX){
            BoxShape box = (BoxShape)b.shape;
            float hx = box.halfExtents.x;
            float hy = box.halfExtents.y;
            float hz = box.halfExtents.z;
            return (float)Math.sqrt(hx*hx + hy*hy + hz*hz);
        }
        AABB a = b.shape.computeAABB(b.position, b.orientation);
        float hx = (a.max.x - a.min.x)*0.5f;
        float hy = (a.max.y - a.min.y)*0.5f;
        float hz = (a.max.z - a.min.z)*0.5f;
        return (float)Math.sqrt(hx*hx + hy*hy + hz*hz);
    }
}
