
package physics3d;

public final class ConvexHullShape implements Shape {
    private final Vec3[] localVerts;
    private final float m;

    public ConvexHullShape(Vec3[] localVerts) {
        if (localVerts == null || localVerts.length < 4)
            throw new IllegalArgumentException("Convex hull needs >= 4 vertices");

        this.localVerts = new Vec3[localVerts.length];
        for (int i = 0; i < localVerts.length; i++)
            this.localVerts[i] = new Vec3(localVerts[i].x, localVerts[i].y, localVerts[i].z);

        this.m = 0.002f;
    }

    @Override public Type type() { return Type.CONVEX_HULL; }

    @Override
    public void supportInto(Vec3 out, Vec3 dirWorld, Vec3 position, Quat orientation) {
        // We assume orientation is (approximately) unit length.
        final float qw = orientation.w, qx = orientation.x, qy = orientation.y, qz = orientation.z;

        // 1) dirLocal = rotateConjugate(q, dirWorld)
        float vx = dirWorld.x, vy = dirWorld.y, vz = dirWorld.z;
        float tx = 2f * (qy * vz - qz * vy);
        float ty = 2f * (qz * vx - qx * vz);
        float tz = 2f * (qx * vy - qy * vx);

        // inverse rotation flips sign on (w*t)
        float lx = vx - qw * tx + (qy * tz - qz * ty);
        float ly = vy - qw * ty + (qz * tx - qx * tz);
        float lz = vz - qw * tz + (qx * ty - qy * tx);

        // 2) Find best vertex in local space along dirLocal
        float best = Float.NEGATIVE_INFINITY;
        float bestX = 0f, bestY = 0f, bestZ = 0f;

        for (Vec3 v : localVerts) {
            float d = v.x * lx + v.y * ly + v.z * lz;
            if (d > best) {
                best = d;
                bestX = v.x; bestY = v.y; bestZ = v.z;
            }
        }

        // 3) Rotate chosen local vertex back to world, add pos
        vx = bestX; vy = bestY; vz = bestZ;
        tx = 2f * (qy * vz - qz * vy);
        ty = 2f * (qz * vx - qx * vz);
        tz = 2f * (qx * vy - qy * vx);

        float wx = vx + qw * tx + (qy * tz - qz * ty);
        float wy = vy + qw * ty + (qz * tx - qx * tz);
        float wz = vz + qw * tz + (qx * ty - qy * tx);

        out.set(position.x + wx, position.y + wy, position.z + wz);
    }

    /** Optional allocating wrapper for compatibility. */
    @Override
    public Vec3 support(Vec3 dirWorld, Vec3 position, Quat orientation) {
        Vec3 out = new Vec3();
        supportInto(out, dirWorld, position, orientation);
        return out;
    }

    /**
     * Allocation-free AABB computation:
     * Writes directly into 'out' (no allocations).
     */
    @Override
    public void computeAABBInto(AABB out, Vec3 position, Quat orientation) {
        out.min.set(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY);
        out.max.set(Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY);

        final float qw = orientation.w, qx = orientation.x, qy = orientation.y, qz = orientation.z;

        for (Vec3 v : localVerts) {
            float vx = v.x, vy = v.y, vz = v.z;

            float tx = 2f * (qy * vz - qz * vy);
            float ty = 2f * (qz * vx - qx * vz);
            float tz = 2f * (qx * vy - qy * vx);

            float rx = vx + qw * tx + (qy * tz - qz * ty);
            float ry = vy + qw * ty + (qz * tx - qx * tz);
            float rz = vz + qw * tz + (qx * ty - qy * tx);

            float wx = position.x + rx;
            float wy = position.y + ry;
            float wz = position.z + rz;

            if (wx < out.min.x) out.min.x = wx;
            if (wy < out.min.y) out.min.y = wy;
            if (wz < out.min.z) out.min.z = wz;

            if (wx > out.max.x) out.max.x = wx;
            if (wy > out.max.y) out.max.y = wy;
            if (wz > out.max.z) out.max.z = wz;
        }

        // Expand by margin without allocating
        out.min.x -= m; out.min.y -= m; out.min.z -= m;
        out.max.x += m; out.max.y += m; out.max.z += m;
    }

    @Override
    public Mat3 inertiaTensor(float mass) {
        // Approximate as local AABB box.
        float minx = Float.POSITIVE_INFINITY, miny = Float.POSITIVE_INFINITY, minz = Float.POSITIVE_INFINITY;
        float maxx = Float.NEGATIVE_INFINITY, maxy = Float.NEGATIVE_INFINITY, maxz = Float.NEGATIVE_INFINITY;

        for (Vec3 v : localVerts) {
            minx = Math.min(minx, v.x); miny = Math.min(miny, v.y); minz = Math.min(minz, v.z);
            maxx = Math.max(maxx, v.x); maxy = Math.max(maxy, v.y); maxz = Math.max(maxz, v.z);
        }

        float hx = (maxx - minx) * 0.5f;
        float hy = (maxy - miny) * 0.5f;
        float hz = (maxz - minz) * 0.5f;

        float Ixx = (mass / 3f) * (hy * hy + hz * hz);
        float Iyy = (mass / 3f) * (hx * hx + hz * hz);
        float Izz = (mass / 3f) * (hx * hx + hy * hy);

        return Mat3.diag(Ixx, Iyy, Izz);
    }

    @Override public float margin() { return m; }
}
