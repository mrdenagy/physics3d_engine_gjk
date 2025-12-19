
package physics3d;

/**
 * Minkowski difference support mapping utilities.
 *
 * IMPORTANT:
 * - supportInto(...) is allocation-free (no new Vec3/SupportPoint created).
 * - support(...) is kept for backwards compatibility but allocates (creates a new SupportPoint).
 */
public final class Minkowski {
    private Minkowski() {}

    public static final class SupportPoint {
        /** Point in Minkowski difference (a - b). */
        public final Vec3 v = new Vec3();
        /** Support point on shape A in world space. */
        public final Vec3 a = new Vec3();
        /** Support point on shape B in world space. */
        public final Vec3 b = new Vec3();

        // Scratch to avoid allocating Vec3.neg(dir)
        private final Vec3 negDir = new Vec3();

        public SupportPoint() {}

        /**
         * Fills this SupportPoint for direction dir (world space).
         * Allocation-free as long as shapes implement Shape.supportInto(...).
         */
        public void set(RigidBody A, RigidBody B, Vec3 dir) {
            // a = support(A, dir)
            A.shape.supportInto(a, dir, A.position, A.orientation);

            // b = support(B, -dir)
            negDir.set(-dir.x, -dir.y, -dir.z);
            B.shape.supportInto(b, negDir, B.position, B.orientation);

            // v = a - b
            v.set(a.x - b.x, a.y - b.y, a.z - b.z);
        }
    }

    /**
     * Allocation-free Minkowski support:
     * Writes the support point into 'out' (no allocations).
     */
    public static void supportInto(RigidBody A, RigidBody B, Vec3 dir, SupportPoint out) {
        out.set(A, B, dir);
    }

    /**
     * Backwards-compatible allocating support:
     * Returns a new SupportPoint each call (allocates).
     *
     * Prefer supportInto(...) in hot loops (GJK/EPA/GJKDistance/CCD).
     */
    public static SupportPoint support(RigidBody A, RigidBody B, Vec3 dir) {
        SupportPoint sp = new SupportPoint();
        sp.set(A, B, dir);
        return sp;
    }
}
