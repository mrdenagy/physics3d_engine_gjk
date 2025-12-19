
package physics3d;

/**
 * Axis-Aligned Bounding Box (AABB).
 * Mutable + allocation-free in use: all operations modify existing instances.
 */
public final class AABB {
    /** Minimum corner. */
    public final Vec3 min = new Vec3();
    /** Maximum corner. */
    public final Vec3 max = new Vec3();

    public AABB() {}

    /** Set from another AABB (no allocations). */
    public AABB set(AABB o) {
        min.set(o.min);
        max.set(o.max);
        return this;
    }

    /** Set from components (no allocations). */
    public AABB set(float minx, float miny, float minz, float maxx, float maxy, float maxz) {
        min.set(minx, miny, minz);
        max.set(maxx, maxy, maxz);
        return this;
    }

    /** Expand by scalar margin (no allocations). */
    public AABB expand(float e) {
        min.x -= e; min.y -= e; min.z -= e;
        max.x += e; max.y += e; max.z += e;
        return this;
    }

    /** Translate by (dx, dy, dz) (no allocations). */
    public AABB translate(float dx, float dy, float dz) {
        min.x += dx; min.y += dy; min.z += dz;
        max.x += dx; max.y += dy; max.z += dz;
        return this;
    }

    /** out = union(a, b) (no allocations). */
    public static void unionInto(AABB out, AABB a, AABB b) {
        out.min.set(
                Math.min(a.min.x, b.min.x),
                Math.min(a.min.y, b.min.y),
                Math.min(a.min.z, b.min.z)
        );
        out.max.set(
                Math.max(a.max.x, b.max.x),
                Math.max(a.max.y, b.max.y),
                Math.max(a.max.z, b.max.z)
        );
    }

    /** True if overlaps other (same semantics as your current method). */
    public boolean overlaps(AABB o) {
        return !(max.x < o.min.x || min.x > o.max.x ||
                 max.y < o.min.y || min.y > o.max.y ||
                 max.z < o.min.z || min.z > o.max.z);
    }

    @Override
    public String toString() {
        return "AABB[min=" + min + ", max=" + max + "]";
    }
}
