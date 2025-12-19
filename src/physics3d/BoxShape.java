
package physics3d;

public final class BoxShape implements Shape {
    public final Vec3 halfExtents;
    private final float m;

    public BoxShape(float hx, float hy, float hz) {
        this.halfExtents = new Vec3(hx, hy, hz);
        this.m = 0.001f;
    }

    @Override public Type type() { return Type.BOX; }

    /**
     * Allocation-free support mapping for an oriented box.
     *
     * Steps:
     * 1) Rotate world direction into LOCAL space using inverse(q) = conjugate(q).
     * 2) Pick extreme corner by sign in local space.
     * 3) Rotate that local point back to WORLD space using q, then add position.
     *
     * No allocations: no new Vec3, no new Quat.
     */
    @Override
    public void supportInto(Vec3 out, Vec3 dirWorld, Vec3 position, Quat orientation) {
        // We assume orientation is (approximately) unit length.
        final float qw = orientation.w, qx = orientation.x, qy = orientation.y, qz = orientation.z;

        // -------- 1) dirLocal = rotateConjugate(q, dirWorld) -----------
        // Standard fast rotate: v' = v + w*t + cross(q.xyz, t), where t = 2*cross(q.xyz, v)
        // For inverse rotation by conjugate(q): v' = v - w*t + cross(q.xyz, t)
        float vx = dirWorld.x, vy = dirWorld.y, vz = dirWorld.z;
        float tx = 2f * (qy * vz - qz * vy);
        float ty = 2f * (qz * vx - qx * vz);
        float tz = 2f * (qx * vy - qy * vx);

        // inverse rotation differs only by the sign on the (w*t) term
        float lx = vx - qw * tx + (qy * tz - qz * ty);
        float ly = vy - qw * ty + (qz * tx - qx * tz);
        float lz = vz - qw * tz + (qx * ty - qy * tx);

        // Pick extreme local corner based on sign of local direction
        float x = (lx >= 0f) ? halfExtents.x : -halfExtents.x;
        float y = (ly >= 0f) ? halfExtents.y : -halfExtents.y;
        float z = (lz >= 0f) ? halfExtents.z : -halfExtents.z;

        // -------- 2) worldPoint = rotate(q, localCorner) + position ----
        vx = x; vy = y; vz = z;
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
     * Allocation-free AABB:
     * Rotates the 8 corners implicitly and tracks min/max.
     * Writes directly into 'out' (no allocations).
     */
    @Override
    public void computeAABBInto(AABB out, Vec3 position, Quat orientation) {
        out.min.set(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY);
        out.max.set(Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY);

        final float qw = orientation.w, qx = orientation.x, qy = orientation.y, qz = orientation.z;

        for (int sx = -1; sx <= 1; sx += 2) {
            for (int sy = -1; sy <= 1; sy += 2) {
                for (int sz = -1; sz <= 1; sz += 2) {
                    float vx = halfExtents.x * sx;
                    float vy = halfExtents.y * sy;
                    float vz = halfExtents.z * sz;

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
            }
        }

        // Expand by margin without allocating
        out.min.x -= m; out.min.y -= m; out.min.z -= m;
        out.max.x += m; out.max.y += m; out.max.z += m;
    }

    @Override
    public Mat3 inertiaTensor(float mass) {
        float hx = halfExtents.x, hy = halfExtents.y, hz = halfExtents.z;

        // solid box about center:
        // Ixx = (1/12)m((2hy)^2 + (2hz)^2) = (m/3)(hy^2 + hz^2)
        float Ixx = (mass / 3f) * (hy * hy + hz * hz);
        float Iyy = (mass / 3f) * (hx * hx + hz * hz);
        float Izz = (mass / 3f) * (hx * hx + hy * hy);

        return Mat3.diag(Ixx, Iyy, Izz);
    }

    @Override public float margin() { return m; }
}
