
package physics3d;

public final class Quat {
    public float w, x, y, z;

    public Quat() { this(1, 0, 0, 0); }
    public Quat(float w, float x, float y, float z) { this.w = w; this.x = x; this.y = y; this.z = z; }

    public static Quat identity() { return new Quat(1, 0, 0, 0); }

    public Quat set(Quat q) { w = q.w; x = q.x; y = q.y; z = q.z; return this; }

    /** Allocating conjugate (kept for compatibility). */
    public Quat conjugate() { return new Quat(w, -x, -y, -z); }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)

    /** Allocation-free conjugate into out. */
    public void conjugateInto(Quat out) { out.w = w; out.x = -x; out.y = -y; out.z = -z; }

    public float dot(Quat q) { return w * q.w + x * q.x + y * q.y + z * q.z; }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)

    public Quat normalized() {
        float n = (float) Math.sqrt(w * w + x * x + y * y + z * z);
        if (n < 1e-8f) return new Quat(1, 0, 0, 0);
        return new Quat(w / n, x / n, y / n, z / n);
    }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)

    /** Allocation-free normalize in place. */
    public Quat normalizeLocal(float eps) {
        float n2 = w*w + x*x + y*y + z*z;
        if (n2 < eps*eps) { w = 1; x = y = z = 0; return this; }
        float inv = 1f / (float)Math.sqrt(n2);
        w *= inv; x *= inv; y *= inv; z *= inv;
        return this;
    }

    public Quat mul(Quat b) {
        return new Quat(
                w * b.w - x * b.x - y * b.y - z * b.z,
                w * b.x + x * b.w + y * b.z - z * b.y,
                w * b.y - x * b.z + y * b.w + z * b.x,
                w * b.z + x * b.y - y * b.x + z * b.w
        );
    }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)

    /** Allocation-free multiply into out: out = this * b */
    public void mulInto(Quat out, Quat b) {
        float ow = w * b.w - x * b.x - y * b.y - z * b.z;
        float ox = w * b.x + x * b.w + y * b.z - z * b.y;
        float oy = w * b.y - x * b.z + y * b.w + z * b.x;
        float oz = w * b.z + x * b.y - y * b.x + z * b.w;
        out.w = ow; out.x = ox; out.y = oy; out.z = oz;
    }

    /**
     * Allocation-free quaternion-vector rotation: out = q * v * q^-1
     * Uses the fast form:
     *   t = 2 * cross(q.xyz, v)
     *   out = v + w*t + cross(q.xyz, t)
     */
    public void rotateInto(Vec3 out, Vec3 v) {
        final float vx = v.x, vy = v.y, vz = v.z;
        final float tx = 2f * (y * vz - z * vy);
        final float ty = 2f * (z * vx - x * vz);
        final float tz = 2f * (x * vy - y * vx);

        out.set(
                vx + w * tx + (y * tz - z * ty),
                vy + w * ty + (z * tx - x * tz),
                vz + w * tz + (x * ty - y * tx)
        );
    }

    /**
     * Allocation-free inverse rotation (conjugate): out = q^-1 * v * q
     * Equivalent to rotating by the conjugate for unit quaternions.
     * Fast form differs only by the sign of (w*t).
     */
    public void rotateConjugateInto(Vec3 out, Vec3 v) {
        final float vx = v.x, vy = v.y, vz = v.z;
        final float tx = 2f * (y * vz - z * vy);
        final float ty = 2f * (z * vx - x * vz);
        final float tz = 2f * (x * vy - y * vx);

        out.set(
                vx - w * tx + (y * tz - z * ty),
                vy - w * ty + (z * tx - x * tz),
                vz - w * tz + (x * ty - y * tx)
        );
    }

    /**
     * Backwards-compatible allocating rotate().
     * Now only allocates the return Vec3 (no temporary Quats). [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)
     */
    public Vec3 rotate(Vec3 v) {
        Vec3 out = new Vec3();
        rotateInto(out, v);
        return out;
    }

    public static Quat fromAxisAngle(Vec3 axis, float angleRad) {
        Vec3 a = axis.normalized();
        float half = angleRad * 0.5f;
        float s = (float) Math.sin(half);
        return new Quat((float) Math.cos(half), a.x * s, a.y * s, a.z * s).normalized();
    }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)

    public static Quat fromAngularVelocity(Vec3 omegaWorld, float dt) {
        float halfDt = 0.5f * dt;
        return new Quat(1f, omegaWorld.x * halfDt, omegaWorld.y * halfDt, omegaWorld.z * halfDt).normalized();
    }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)

    /**
     * Returns the relative rotation qRel = qB * inverse(qA) (assuming unit quaternions).
     */
    public static Quat relative(Quat qA, Quat qB) {
        return qB.mul(qA.conjugate()).normalized();
    }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)

    /**
     * Extract twist angle around a given world axis from a relative quaternion.
     * Uses swing-twist decomposition by projecting vector part onto the axis.
     * The returned angle is in [-pi, pi].
     */
    public static float twistAngleAroundAxis(Quat qRel, Vec3 axisWorld) {
        Vec3 a = axisWorld.normalized();
        float proj = qRel.x * a.x + qRel.y * a.y + qRel.z * a.z;
        Quat twist = new Quat(qRel.w, a.x * proj, a.y * proj, a.z * proj).normalized();

        float sinHalf = (float) Math.sqrt(twist.x * twist.x + twist.y * twist.y + twist.z * twist.z);
        float angle = 2f * (float) Math.atan2(sinHalf, twist.w);

        float sign = (proj >= 0f) ? 1f : -1f;
        angle *= sign;

        while (angle > Math.PI) angle -= (float) (2 * Math.PI);
        while (angle < -Math.PI) angle += (float) (2 * Math.PI);

        return angle;
    }  // [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Quat.java)
}
