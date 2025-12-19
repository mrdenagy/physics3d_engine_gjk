
package physics3d;

public interface Shape {
    enum Type { SPHERE, BOX, CONVEX_HULL }

    Type type();

    /**
     * Allocation-free support mapping.
     * Writes the farthest point in the direction (dirWorld) into 'out'.
     */
    void supportInto(Vec3 out, Vec3 dirWorld, Vec3 position, Quat orientation);

    /**
     * Backwards-compatible allocating wrapper.
     * Prefer supportInto(...) in hot loops (GJK/CCD) to avoid allocations.
     */
    default Vec3 support(Vec3 dirWorld, Vec3 position, Quat orientation) {
        Vec3 out = new Vec3();
        supportInto(out, dirWorld, position, orientation);
        return out;
    }

    /**
     * Allocation-free AABB computation in world space.
     * Writes the AABB into 'out' (no allocations).
     */
    void computeAABBInto(AABB out, Vec3 position, Quat orientation);

    /**
     * Backwards-compatible allocating wrapper.
     * Prefer computeAABBInto(...) in broadphase/hot loops to avoid allocations.
     */
    default AABB computeAABB(Vec3 position, Quat orientation) {
        AABB a = new AABB();
        computeAABBInto(a, position, orientation);
        return a;
    }

    Mat3 inertiaTensor(float mass);

    /** Small padding for robustness (like a "skin" radius). */
    float margin();
}
