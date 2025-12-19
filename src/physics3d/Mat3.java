
package physics3d;

/**
 * 3x3 matrix (row-major).
 *
 * Existing allocating API preserved.
 * Added allocation-free "Into" methods for hot loops.
 */
public final class Mat3 {
    public float m00, m01, m02,
                 m10, m11, m12,
                 m20, m21, m22;

    public Mat3() { setIdentity(); }

    public Mat3(float m00, float m01, float m02,
                float m10, float m11, float m12,
                float m20, float m21, float m22) {
        this.m00 = m00; this.m01 = m01; this.m02 = m02;
        this.m10 = m10; this.m11 = m11; this.m12 = m12;
        this.m20 = m20; this.m21 = m21; this.m22 = m22;
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public Mat3 setIdentity() {
        m00 = 1; m01 = 0; m02 = 0;
        m10 = 0; m11 = 1; m12 = 0;
        m20 = 0; m21 = 0; m22 = 1;
        return this;
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public Mat3 set(Mat3 o) {
        m00=o.m00; m01=o.m01; m02=o.m02;
        m10=o.m10; m11=o.m11; m12=o.m12;
        m20=o.m20; m21=o.m21; m22=o.m22;
        return this;
    }

    public static Mat3 diag(float a, float b, float c) {
        return new Mat3(a, 0, 0,
                        0, b, 0,
                        0, 0, c);
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    // -------------------------
    // Allocating methods (kept)
    // -------------------------
    public Vec3 mul(Vec3 v) {
        return new Vec3(
                m00 * v.x + m01 * v.y + m02 * v.z,
                m10 * v.x + m11 * v.y + m12 * v.z,
                m20 * v.x + m21 * v.y + m22 * v.z
        );
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public Mat3 mul(Mat3 b) {
        Mat3 r = new Mat3();
        r.m00 = m00 * b.m00 + m01 * b.m10 + m02 * b.m20;
        r.m01 = m00 * b.m01 + m01 * b.m11 + m02 * b.m21;
        r.m02 = m00 * b.m02 + m01 * b.m12 + m02 * b.m22;

        r.m10 = m10 * b.m00 + m11 * b.m10 + m12 * b.m20;
        r.m11 = m10 * b.m01 + m11 * b.m11 + m12 * b.m21;
        r.m12 = m10 * b.m02 + m11 * b.m12 + m12 * b.m22;

        r.m20 = m20 * b.m00 + m21 * b.m10 + m22 * b.m20;
        r.m21 = m20 * b.m01 + m21 * b.m11 + m22 * b.m21;
        r.m22 = m20 * b.m02 + m21 * b.m12 + m22 * b.m22;
        return r;
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public Mat3 add(Mat3 b) {
        return new Mat3(
                m00 + b.m00, m01 + b.m01, m02 + b.m02,
                m10 + b.m10, m11 + b.m11, m12 + b.m12,
                m20 + b.m20, m21 + b.m21, m22 + b.m22
        );
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public Mat3 sub(Mat3 b) {
        return new Mat3(
                m00 - b.m00, m01 - b.m01, m02 - b.m02,
                m10 - b.m10, m11 - b.m11, m12 - b.m12,
                m20 - b.m20, m21 - b.m21, m22 - b.m22
        );
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public Mat3 transpose() {
        return new Mat3(
                m00, m10, m20,
                m01, m11, m21,
                m02, m12, m22
        );
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public float det() {
        return m00 * (m11 * m22 - m12 * m21)
             - m01 * (m10 * m22 - m12 * m20)
             + m02 * (m10 * m21 - m11 * m20);
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public Mat3 inverse() {
        float d = det();
        if (Math.abs(d) < 1e-12f) return new Mat3().setIdentity();
        float inv = 1f / d;

        float i00 = (m11 * m22 - m12 * m21) * inv;
        float i01 = (m02 * m21 - m01 * m22) * inv;
        float i02 = (m01 * m12 - m02 * m11) * inv;

        float i10 = (m12 * m20 - m10 * m22) * inv;
        float i11 = (m00 * m22 - m02 * m20) * inv;
        float i12 = (m02 * m10 - m00 * m12) * inv;

        float i20 = (m10 * m21 - m11 * m20) * inv;
        float i21 = (m01 * m20 - m00 * m21) * inv;
        float i22 = (m00 * m11 - m01 * m10) * inv;

        return new Mat3(i00, i01, i02,
                        i10, i11, i12,
                        i20, i21, i22);
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public static Mat3 fromQuat(Quat q) {
        float w = q.w, x = q.x, y = q.y, z = q.z;
        float xx = x * x, yy = y * y, zz = z * z;
        float xy = x * y, xz = x * z, yz = y * z;
        float wx = w * x, wy = w * y, wz = w * z;

        return new Mat3(
                1 - 2 * (yy + zz), 2 * (xy - wz),     2 * (xz + wy),
                2 * (xy + wz),     1 - 2 * (xx + zz), 2 * (yz - wx),
                2 * (xz - wy),     2 * (yz + wx),     1 - 2 * (xx + yy)
        );
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    /** Skew-symmetric cross-product matrix [v]_x such that [v]_x * w = v x w */
    public static Mat3 skew(Vec3 v) {
        return new Mat3(
                0, -v.z, v.y,
                v.z, 0, -v.x,
                -v.y, v.x, 0
        );
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    // -----------------------------------
    // Allocation-free "Into" additions
    // -----------------------------------
    public void mulInto(Vec3 out, Vec3 v) {
        out.set(
                m00 * v.x + m01 * v.y + m02 * v.z,
                m10 * v.x + m11 * v.y + m12 * v.z,
                m20 * v.x + m21 * v.y + m22 * v.z
        );
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public void addInto(Mat3 out, Mat3 b) {
        out.m00 = m00 + b.m00; out.m01 = m01 + b.m01; out.m02 = m02 + b.m02;
        out.m10 = m10 + b.m10; out.m11 = m11 + b.m11; out.m12 = m12 + b.m12;
        out.m20 = m20 + b.m20; out.m21 = m21 + b.m21; out.m22 = m22 + b.m22;
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public void subInto(Mat3 out, Mat3 b) {
        out.m00 = m00 - b.m00; out.m01 = m01 - b.m01; out.m02 = m02 - b.m02;
        out.m10 = m10 - b.m10; out.m11 = m11 - b.m11; out.m12 = m12 - b.m12;
        out.m20 = m20 - b.m20; out.m21 = m21 - b.m21; out.m22 = m22 - b.m22;
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public void transposeInto(Mat3 out) {
        out.m00 = m00; out.m01 = m10; out.m02 = m20;
        out.m10 = m01; out.m11 = m11; out.m12 = m21;
        out.m20 = m02; out.m21 = m12; out.m22 = m22;
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)

    public void mulInto(Mat3 out, Mat3 b) {
        float r00 = m00 * b.m00 + m01 * b.m10 + m02 * b.m20;
        float r01 = m00 * b.m01 + m01 * b.m11 + m02 * b.m21;
        float r02 = m00 * b.m02 + m01 * b.m12 + m02 * b.m22;

        float r10 = m10 * b.m00 + m11 * b.m10 + m12 * b.m20;
        float r11 = m10 * b.m01 + m11 * b.m11 + m12 * b.m21;
        float r12 = m10 * b.m02 + m11 * b.m12 + m12 * b.m22;

        float r20 = m20 * b.m00 + m21 * b.m10 + m22 * b.m20;
        float r21 = m20 * b.m01 + m21 * b.m11 + m22 * b.m21;
        float r22 = m20 * b.m02 + m21 * b.m12 + m22 * b.m22;

        out.m00 = r00; out.m01 = r01; out.m02 = r02;
        out.m10 = r10; out.m11 = r11; out.m12 = r12;
        out.m20 = r20; out.m21 = r21; out.m22 = r22;
    }  // [2](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/Mat3.java)
}
