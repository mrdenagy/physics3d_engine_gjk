
package physics3d;

public final class RigidBody {
    public enum Type { STATIC, DYNAMIC }

    public final Shape shape;
    public Type type;

    public Vec3 position = new Vec3();
    public Quat orientation = Quat.identity();

    public Vec3 linearVelocity = new Vec3();
    public Vec3 angularVelocity = new Vec3(); // world-space omega

    public Vec3 force = new Vec3();
    public Vec3 torque = new Vec3();

    public float restitution = 0.2f;
    public float friction = 0.6f;

    // Damping (per-second). Typical: 0.0 to 5.0
    public float linearDamping = 0.0f;
    public float angularDamping = 0.0f;

    public float mass, invMass;
    public Mat3 inertiaBody, invInertiaBody; // invInertiaBody is diagonal in this engine

    // Sleeping
    public boolean canSleep = true;
    public boolean awake = true;
    public float sleepTime = 0f;

    // Continuous collision detection (treat as fast-moving bullet)
    public boolean bullet = false;

    // ✅ Cached AABB (allocation-free broadphase integration)
    public final AABB cachedAabb = new AABB();

    // ---- Scratch vectors (avoid per-call allocations) ----
    private final Vec3 tmpR = new Vec3();       // worldPoint - position
    private final Vec3 tmpCross = new Vec3();   // cross products
    private final Vec3 tmpDw = new Vec3();      // angular delta
    private final Vec3 tmpRot = new Vec3();     // rotated vector
    private final Vec3 tmpBody = new Vec3();    // body-space vector
    private final Vec3 tmpScaled = new Vec3();  // scaled body-space vector

    public RigidBody(Shape shape, Type type, float mass, Vec3 position) {
        this.shape = shape;
        this.type = type;
        this.position.set(position);

        if (type == Type.STATIC || mass <= 0f) {
            this.mass = Float.POSITIVE_INFINITY;
            this.invMass = 0f;
            this.inertiaBody = Mat3.diag(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY);
            this.invInertiaBody = Mat3.diag(0, 0, 0);
            this.awake = false;
            this.canSleep = false;
        } else {
            this.mass = mass;
            this.invMass = 1f / mass;
            this.inertiaBody = shape.inertiaTensor(mass);
            this.invInertiaBody = Mat3.diag(
                    MathUtil.safeInv(inertiaBody.m00),
                    MathUtil.safeInv(inertiaBody.m11),
                    MathUtil.safeInv(inertiaBody.m22)
            );
            this.awake = true;
            this.canSleep = true;
        }

        // ✅ Initialize cached AABB once (allocation-free)
        updateAABB();
    }

    public boolean isStatic() {
        return invMass == 0f || type == Type.STATIC;
    }

    public void setAwake(boolean awake) {
        if (isStatic()) {
            this.awake = false;
            return;
        }

        if (awake) {
            this.awake = true;
            this.sleepTime = 0f;
        } else {
            this.awake = false;
            this.sleepTime = 0f;
            this.linearVelocity.set(0, 0, 0);
            this.angularVelocity.set(0, 0, 0);
            this.force.set(0, 0, 0);
            this.torque.set(0, 0, 0);
        }
    }

    public void applyForce(Vec3 f) {
        if (isStatic()) return;
        if (!awake) setAwake(true);
        force.addLocal(f);
    }

    /**
     * Allocation-free: uses scratch vectors instead of worldPoint.sub(position) and r.cross(f).
     */
    public void applyForceAtPoint(Vec3 f, Vec3 worldPoint) {
        if (isStatic()) return;
        if (!awake) setAwake(true);

        force.addLocal(f);

        // r = worldPoint - position (no allocation)
        tmpR.set(worldPoint.x - position.x, worldPoint.y - position.y, worldPoint.z - position.z);

        // torque += r x f (no allocation)
        tmpCross.setCross(tmpR, f);
        torque.addLocal(tmpCross);
    }

    /**
     * Allocation-free impulse application:
     * - linearVelocity += impulse * invMass
     * - angularVelocity += invI_world * (r x impulse)
     *
     * Avoids allocating Vec3 from mul/sub/cross and avoids allocating Mat3 by doing:
     * invI_world * v = R( invI_body * (R^T v) ), with invI_body diagonal.
     */
    public void applyImpulse(Vec3 impulse, Vec3 worldPoint) {
        if (isStatic()) return;
        if (!awake) setAwake(true);

        // linearVelocity += impulse * invMass (no allocation)
        linearVelocity.addScaledLocal(impulse, invMass);

        // r = worldPoint - position
        tmpR.set(worldPoint.x - position.x, worldPoint.y - position.y, worldPoint.z - position.z);

        // tmpCross = r x impulse
        tmpCross.setCross(tmpR, impulse);

        // tmpDw = invI_world * tmpCross (no allocation)
        invInertiaWorldMulInto(tmpDw, tmpCross);

        angularVelocity.addLocal(tmpDw);
    }

    /**
     * Original method kept (may allocate Mat3 depending on your Mat3 implementation).
     * Prefer invInertiaWorldMulInto(...) in hot paths.
     */
    public Mat3 invInertiaWorld() {
        Mat3 R = Mat3.fromQuat(orientation);
        return R.mul(invInertiaBody).mul(R.transpose());
    }

    // ---------------------------------------------------------------------
    // ✅ Cached AABB API (near-zero GC)
    // ---------------------------------------------------------------------

    /**
     * Updates cachedAabb in-place (allocation-free).
     * Call this after changing position/orientation.
     */
    public void updateAABB() {
        shape.computeAABBInto(cachedAabb, position, orientation);
    }

    /**
     * Returns a reference to the cached AABB (no allocations).
     * This replaces the old allocating behavior.
     */
    public AABB aabb() {
        return cachedAabb;
    }

    /**
     * Backwards-compatible allocating AABB if needed in non-hot code.
     * Prefer aabb() + updateAABB() for near-zero GC.
     */
    public AABB aabbAllocating() {
        return shape.computeAABB(position, orientation);
    }

    // ---------------------------------------------------------------------
    // Allocation-free coordinate transforms (avoid Quat.rotate allocations)
    // ---------------------------------------------------------------------

    /** Allocation-free: out = position + rotate(orientation, local) */
    public void worldPointFromLocalInto(Vec3 out, Vec3 local) {
        rotateInto(tmpRot, orientation, local); // tmpRot = R * local
        out.set(position.x + tmpRot.x, position.y + tmpRot.y, position.z + tmpRot.z);
    }

    /** Backwards-compatible allocating wrapper (kept). */
    public Vec3 worldPointFromLocal(Vec3 local) {
        Vec3 out = new Vec3();
        worldPointFromLocalInto(out, local);
        return out;
    }

    /** Allocation-free: out = rotateConjugate(orientation, world - position) */
    public void localPointFromWorldInto(Vec3 out, Vec3 world) {
        tmpR.set(world.x - position.x, world.y - position.y, world.z - position.z);
        rotateConjugateInto(out, orientation, tmpR);
    }

    /** Backwards-compatible allocating wrapper (kept). */
    public Vec3 localPointFromWorld(Vec3 world) {
        Vec3 out = new Vec3();
        localPointFromWorldInto(out, world);
        return out;
    }

    // ---------------------------------------------------------------------
    // Allocation-free inv inertia world multiply (no Mat3 allocations)
    // invI_world * v = R( invI_body * (R^T v) ), invI_body diagonal
    // ---------------------------------------------------------------------
    private void invInertiaWorldMulInto(Vec3 out, Vec3 vWorld) {
        // vBody = R^T * vWorld (inverse rotation)
        rotateConjugateInto(tmpBody, orientation, vWorld);

        // scale by diagonal inv inertia in body space
        tmpScaled.set(
                tmpBody.x * invInertiaBody.m00,
                tmpBody.y * invInertiaBody.m11,
                tmpBody.z * invInertiaBody.m22
        );

        // out = R * tmpScaled (forward rotation)
        rotateInto(out, orientation, tmpScaled);
    }

    // ---------------------------------------------------------------------
    // Quaternion rotate helpers (allocation-free, float math)
    // ---------------------------------------------------------------------
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

        // Inverse rotation (conjugate): flip sign on qw*t term
        out.set(
                vx - qw * tx + (qy * tz - qz * ty),
                vy - qw * ty + (qz * tx - qx * tz),
                vz - qw * tz + (qx * ty - qy * tx)
        );
    }
}
