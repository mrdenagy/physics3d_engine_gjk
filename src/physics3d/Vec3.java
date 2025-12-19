
package physics3d;

/**
 * 3D vector with both allocating convenience methods (add/sub/mul/cross)
 * and allocation-free "Into"/"set*" utilities for hot loops.
 *
 * Existing API preserved; additional methods added for near-zero GC physics loops.
 */
public final class Vec3 {
    public float x, y, z;

    public Vec3() { this(0, 0, 0); }
    public Vec3(float x, float y, float z) { this.x = x; this.y = y; this.z = z; }

    public Vec3 set(float x, float y, float z) { this.x = x; this.y = y; this.z = z; return this; }
    public Vec3 set(Vec3 v) { this.x = v.x; this.y = v.y; this.z = v.z; return this; }

    // ------------------------------------------------------------------
    // Allocating convenience methods (kept for compatibility / readability)
    // ------------------------------------------------------------------
    public Vec3 add(Vec3 v) { return new Vec3(x + v.x, y + v.y, z + v.z); }
    public Vec3 sub(Vec3 v) { return new Vec3(x - v.x, y - v.y, z - v.z); }
    public Vec3 mul(float s) { return new Vec3(x * s, y * s, z * s); }

    public Vec3 cross(Vec3 v) {
        return new Vec3(
                y * v.z - z * v.y,
                z * v.x - x * v.z,
                x * v.y - y * v.x
        );
    }

    public Vec3 abs() { return new Vec3(Math.abs(x), Math.abs(y), Math.abs(z)); }
    public static Vec3 neg(Vec3 v) { return new Vec3(-v.x, -v.y, -v.z); }

    public Vec3 normalized() {
        Vec3 out = new Vec3();
        normalizeInto(out, this, 1e-8f);
        return out;
    }

    // ------------------------------------------------------------------
    // In-place local mutation (allocation-free)
    // ------------------------------------------------------------------
    public void addLocal(Vec3 v) { x += v.x; y += v.y; z += v.z; }
    public void subLocal(Vec3 v) { x -= v.x; y -= v.y; z -= v.z; }
    public void mulLocal(float s) { x *= s; y *= s; z *= s; }

    /** this += v * s (very common in solvers / integration). */
    public void addScaledLocal(Vec3 v, float s) { x += v.x * s; y += v.y * s; z += v.z * s; }

    /** this -= v * s */
    public void subScaledLocal(Vec3 v, float s) { x -= v.x * s; y -= v.y * s; z -= v.z * s; }

    /** this = -this */
    public void negLocal() { x = -x; y = -y; z = -z; }

    /** this = abs(this) */
    public void absLocal() { x = Math.abs(x); y = Math.abs(y); z = Math.abs(z); }

    /** Normalize this in place with epsilon threshold. */
    public Vec3 normalizeLocal(float eps) {
        float l2 = x * x + y * y + z * z;
        if (l2 < eps * eps) { x = y = z = 0f; return this; }
        float inv = 1f / (float) Math.sqrt(l2);
        x *= inv; y *= inv; z *= inv;
        return this;
    }

    // ------------------------------------------------------------------
    // Scalar ops
    // ------------------------------------------------------------------
    public float dot(Vec3 v) { return x * v.x + y * v.y + z * v.z; }
    public static float dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

    public float len2() { return x * x + y * y + z * z; }
    public float len() { return (float) Math.sqrt(len2()); }

    // ------------------------------------------------------------------
    // Allocation-free "Into" utilities (static)
    // (You already had many; kept and extended.)
    // ------------------------------------------------------------------
    public static void addInto(Vec3 out, Vec3 a, Vec3 b) {
        out.x = a.x + b.x; out.y = a.y + b.y; out.z = a.z + b.z;
    }

    public static void subInto(Vec3 out, Vec3 a, Vec3 b) {
        out.x = a.x - b.x; out.y = a.y - b.y; out.z = a.z - b.z;
    }

    public static void mulInto(Vec3 out, Vec3 a, float s) {
        out.x = a.x * s; out.y = a.y * s; out.z = a.z * s;
    }

    /** out = a + b*s */
    public static void maddInto(Vec3 out, Vec3 a, Vec3 b, float s) {
        out.x = a.x + b.x * s;
        out.y = a.y + b.y * s;
        out.z = a.z + b.z * s;
    }

    /** out = a - b*s */
    public static void msubInto(Vec3 out, Vec3 a, Vec3 b, float s) {
        out.x = a.x - b.x * s;
        out.y = a.y - b.y * s;
        out.z = a.z - b.z * s;
    }

    public static void negInto(Vec3 out, Vec3 v) {
        out.x = -v.x; out.y = -v.y; out.z = -v.z;
    }

    public static void absInto(Vec3 out, Vec3 v) {
        out.x = Math.abs(v.x); out.y = Math.abs(v.y); out.z = Math.abs(v.z);
    }

    public static void crossInto(Vec3 out, Vec3 a, Vec3 b) {
        float cx = a.y * b.z - a.z * b.y;
        float cy = a.z * b.x - a.x * b.z;
        float cz = a.x * b.y - a.y * b.x;
        out.x = cx; out.y = cy; out.z = cz;
    }

    /** out = normalize(v) with epsilon. */
    public static void normalizeInto(Vec3 out, Vec3 v, float eps) {
        float l2 = v.x * v.x + v.y * v.y + v.z * v.z;
        if (l2 < eps * eps) { out.x = out.y = out.z = 0f; return; }
        float inv = 1f / (float) Math.sqrt(l2);
        out.x = v.x * inv; out.y = v.y * inv; out.z = v.z * inv;
    }

    // ------------------------------------------------------------------
    // Allocation-free "set*" instance utilities (new, very useful)
    // ------------------------------------------------------------------
    /** this = a + b */
    public Vec3 setAdd(Vec3 a, Vec3 b) { x = a.x + b.x; y = a.y + b.y; z = a.z + b.z; return this; }

    /** this = a - b */
    public Vec3 setSub(Vec3 a, Vec3 b) { x = a.x - b.x; y = a.y - b.y; z = a.z - b.z; return this; }

    /** this = a * s */
    public Vec3 setMul(Vec3 a, float s) { x = a.x * s; y = a.y * s; z = a.z * s; return this; }

    /** this = -a */
    public Vec3 setNeg(Vec3 a) { x = -a.x; y = -a.y; z = -a.z; return this; }

    /** this = abs(a) */
    public Vec3 setAbs(Vec3 a) { x = Math.abs(a.x); y = Math.abs(a.y); z = Math.abs(a.z); return this; }

    /** this = a x b */
    public Vec3 setCross(Vec3 a, Vec3 b) {
        float cx = a.y * b.z - a.z * b.y;
        float cy = a.z * b.x - a.x * b.z;
        float cz = a.x * b.y - a.y * b.x;
        x = cx; y = cy; z = cz;
        return this;
    }

    /** this = a + b*s */
    public Vec3 setMadd(Vec3 a, Vec3 b, float s) {
        x = a.x + b.x * s;
        y = a.y + b.y * s;
        z = a.z + b.z * s;
        return this;
    }

    /** this = normalize(v) */
    public Vec3 setNormalized(Vec3 v, float eps) {
        float l2 = v.x * v.x + v.y * v.y + v.z * v.z;
        if (l2 < eps * eps) { x = y = z = 0f; return this; }
        float inv = 1f / (float) Math.sqrt(l2);
        x = v.x * inv; y = v.y * inv; z = v.z * inv;
        return this;
    }

    @Override
    public String toString() {
        return "Vec3(" + x + "," + y + "," + z + ")";
    }
}