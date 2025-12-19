
package physics3d;

public final class SphereShape implements Shape {
    public final float r;

    public SphereShape(float r) { this.r = r; }

    @Override public Type type() { return Type.SPHERE; }

    /**
     * Allocation-free support mapping for a sphere.
     * out = position + normalize(dirWorld) * r
     */
    @Override
    public void supportInto(Vec3 out, Vec3 dirWorld, Vec3 position, Quat orientation) {
        float x = dirWorld.x, y = dirWorld.y, z = dirWorld.z;
        float len2 = x * x + y * y + z * z;

        // If dir is near zero, choose an arbitrary direction (+X).
        if (len2 < 1e-12f) {
            out.set(position.x + r, position.y, position.z);
            return;
        }

        float invLen = 1f / (float) Math.sqrt(len2);
        out.set(
                position.x + x * invLen * r,
                position.y + y * invLen * r,
                position.z + z * invLen * r
        );
    }

    /**
     * Optional allocating wrapper for compatibility.
     * Prefer supportInto(...) in hot loops to avoid allocations. [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/SphereShape.java)
     */
    @Override
    public Vec3 support(Vec3 dirWorld, Vec3 position, Quat orientation) {
        Vec3 out = new Vec3();
        supportInto(out, dirWorld, position, orientation);
        return out;
    }

    /**
     * Allocation-free AABB computation for a sphere (orientation irrelevant).
     * Writes directly into 'out' (no allocations). [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/SphereShape.java)
     */
    @Override
    public void computeAABBInto(AABB out, Vec3 position, Quat orientation) {
        out.min.set(position.x - r, position.y - r, position.z - r);
        out.max.set(position.x + r, position.y + r, position.z + r);
    }

    @Override
    public Mat3 inertiaTensor(float mass) {
        // solid sphere: I = 2/5 m r^2
        float I = 0.4f * mass * r * r;
        return Mat3.diag(I, I, I);
    }

    @Override public float margin() { return 0.001f; }
}
