
package physics3d;

import java.util.ArrayList;

/**
 * GJK intersection test for convex shapes in 3D.
 *
 * Optimized:
 * - Allocation-free in the hot loop (fixed simplex, reused SupportPoints/Vec3).
 * - Uses Minkowski.supportInto(...) to avoid support-mapping allocations.
 *
 * NOTE:
 * Result.simplex is still returned as an ArrayList for compatibility with EPA.expand(...).
 * That list is populated from the fixed simplex at the end (small, bounded work).
 */
public final class GJK {
    private GJK() {}

    public static final class Result {
        public boolean intersect;
        public ArrayList<Minkowski.SupportPoint> simplex = new ArrayList<>(4);
        public Vec3 direction = new Vec3(1, 0, 0);
    }

    // Workspace for allocation-free hot loop.
    private static final class Workspace {
        // Fixed simplex (newest is index 0)
        final Minkowski.SupportPoint[] s = new Minkowski.SupportPoint[] {
                new Minkowski.SupportPoint(),
                new Minkowski.SupportPoint(),
                new Minkowski.SupportPoint(),
                new Minkowski.SupportPoint()
        };
        int n = 0;

        final Vec3 d = new Vec3(1, 0, 0);
        final Vec3 AO = new Vec3();
        final Vec3 AB = new Vec3();
        final Vec3 AC = new Vec3();
        final Vec3 AD = new Vec3();

        final Vec3 ABC = new Vec3();
        final Vec3 ACD = new Vec3();
        final Vec3 ADB = new Vec3();

        final Vec3 tmp1 = new Vec3();
        final Vec3 tmp2 = new Vec3();
        final Vec3 tmp3 = new Vec3();
    }

    private static final ThreadLocal<Workspace> TLS = ThreadLocal.withInitial(Workspace::new);

    public static Result intersect(RigidBody A, RigidBody B) {
        Workspace ws = TLS.get();
        ws.n = 0;

        // Initial direction from A to B
        ws.d.set(B.position.x - A.position.x, B.position.y - A.position.y, B.position.z - A.position.z);
        if (ws.d.len2() < 1e-8f) ws.d.set(1, 0, 0);

        // First support point
        Minkowski.supportInto(A, B, ws.d, ws.s[0]);
        ws.n = 1;

        // New direction toward origin: d = -v
        ws.d.set(-ws.s[0].v.x, -ws.s[0].v.y, -ws.s[0].v.z);

        Result r = new Result();

        for (int iter = 0; iter < 30; iter++) {
            // Get next support point in direction d
            // Insert new point at index 0: shift simplex right
            shiftRight(ws);
            Minkowski.supportInto(A, B, ws.d, ws.s[0]);
            if (ws.n < 4) ws.n++;

            // If the newest point isn't past the origin in direction d, no intersection
            if (ws.s[0].v.dot(ws.d) < 0f) {
                r.intersect = false;
                r.direction.set(ws.d);
                r.simplex = exportSimplex(ws);
                return r;
            }

            // Update simplex and direction
            if (containsOrigin(ws)) {
                r.intersect = true;
                r.direction.set(ws.d);
                r.simplex = exportSimplex(ws);
                return r;
            }
        }

        r.intersect = false;
        r.direction.set(ws.d);
        r.simplex = exportSimplex(ws);
        return r;
    }

    /** Export fixed simplex into ArrayList for EPA compatibility (bounded small allocation). */
    private static ArrayList<Minkowski.SupportPoint> exportSimplex(Workspace ws) {
        ArrayList<Minkowski.SupportPoint> out = new ArrayList<>(ws.n);
        for (int i = 0; i < ws.n; i++) out.add(ws.s[i]);
        return out;
    }

    /** Shift simplex entries right by one (keeping max size 4). */
    private static void shiftRight(Workspace ws) {
        // Move 2->3, 1->2, 0->1 (drop 3)
        ws.s[3] = ws.s[2];
        ws.s[2] = ws.s[1];
        ws.s[1] = ws.s[0];
        // ws.s[0] will be overwritten by supportInto; but we must keep an object there.
        // Reuse the last one we dropped to keep four persistent objects:
        ws.s[0] = (ws.s[3] != null) ? ws.s[3] : new Minkowski.SupportPoint();
    }

    /**
     * Updates ws.n and ws.d based on current simplex.
     * Returns true if origin is inside the simplex (intersection).
     */
    private static boolean containsOrigin(Workspace ws) {
        // A is newest point at s[0]
        Minkowski.SupportPoint A = ws.s[0];
        ws.AO.set(-A.v.x, -A.v.y, -A.v.z);

        if (ws.n == 1) {
            ws.d.set(ws.AO);
            return false;
        }

        if (ws.n == 2) {
            Minkowski.SupportPoint B = ws.s[1];

            // AB = B - A
            ws.AB.set(B.v.x - A.v.x, B.v.y - A.v.y, B.v.z - A.v.z);

            // d = tripleCross(AB, AO, AB) (perpendicular toward origin)
            tripleCrossInto(ws.AB, ws.AO, ws.AB, ws.tmp1);
            if (ws.tmp1.len2() < 1e-10f) {
                perpendicularInto(ws.AB, ws.tmp1);
            }
            ws.d.set(ws.tmp1);
            return false;
        }

        if (ws.n == 3) {
            Minkowski.SupportPoint B = ws.s[1];
            Minkowski.SupportPoint C = ws.s[2];

            // AB, AC
            ws.AB.set(B.v.x - A.v.x, B.v.y - A.v.y, B.v.z - A.v.z);
            ws.AC.set(C.v.x - A.v.x, C.v.y - A.v.y, C.v.z - A.v.z);

            // ABC = AB x AC
            crossInto(ws.AB, ws.AC, ws.ABC);

            // abPerp = tripleCross(AC, AB, AB)
            tripleCrossInto(ws.AC, ws.AB, ws.AB, ws.tmp1);
            if (ws.tmp1.dot(ws.AO) > 0f) {
                // Remove C -> keep (A,B)
                ws.s[2] = ws.s[1];
                ws.s[1] = ws.s[0];
                ws.n = 2;
                ws.d.set(ws.tmp1);
                return false;
            }

            // acPerp = tripleCross(AB, AC, AC)
            tripleCrossInto(ws.AB, ws.AC, ws.AC, ws.tmp2);
            if (ws.tmp2.dot(ws.AO) > 0f) {
                // Remove B -> keep (A,C)
                ws.s[1] = ws.s[0];
                ws.n = 2;
                ws.d.set(ws.tmp2);
                return false;
            }

            // Origin is above or below triangle
            if (ws.ABC.dot(ws.AO) > 0f) {
                ws.d.set(ws.ABC);
            } else {
                // Swap B and C (keep winding) and flip normal
                Minkowski.SupportPoint tmp = ws.s[1];
                ws.s[1] = ws.s[2];
                ws.s[2] = tmp;
                ws.d.set(-ws.ABC.x, -ws.ABC.y, -ws.ABC.z);
            }
            return false;
        }

        // ws.n == 4 (tetrahedron)
        Minkowski.SupportPoint B = ws.s[1];
        Minkowski.SupportPoint C = ws.s[2];
        Minkowski.SupportPoint D = ws.s[3];

        // AB, AC, AD
        ws.AB.set(B.v.x - A.v.x, B.v.y - A.v.y, B.v.z - A.v.z);
        ws.AC.set(C.v.x - A.v.x, C.v.y - A.v.y, C.v.z - A.v.z);
        ws.AD.set(D.v.x - A.v.x, D.v.y - A.v.y, D.v.z - A.v.z);

        // Face normals
        crossInto(ws.AB, ws.AC, ws.ABC);
        crossInto(ws.AC, ws.AD, ws.ACD);
        crossInto(ws.AD, ws.AB, ws.ADB);

        // If origin is outside any face, drop opposite point and continue toward that face
        if (ws.ABC.dot(ws.AO) > 0f) {
            // Remove D -> keep (A,B,C)
            ws.n = 3;
            ws.d.set(ws.ABC);
            return false;
        }
        if (ws.ACD.dot(ws.AO) > 0f) {
            // Remove B -> keep (A,C,D) : reorder to (A,C,D)
            ws.s[1] = ws.s[2];
            ws.s[2] = ws.s[3];
            ws.n = 3;
            ws.d.set(ws.ACD);
            return false;
        }
        if (ws.ADB.dot(ws.AO) > 0f) {
            // Remove C -> keep (A,D,B) : reorder to (A,D,B)
            Minkowski.SupportPoint oldB = ws.s[1];
            ws.s[1] = ws.s[3];
            ws.s[2] = oldB;
            ws.n = 3;
            ws.d.set(ws.ADB);
            return false;
        }

        // Origin inside tetrahedron
        return true;
    }

    /** out = a x b */
    private static void crossInto(Vec3 a, Vec3 b, Vec3 out) {
        float x = a.y * b.z - a.z * b.y;
        float y = a.z * b.x - a.x * b.z;
        float z = a.x * b.y - a.y * b.x;
        out.set(x, y, z);
    }

    /** out = a x (b x c) */
    private static void tripleCrossInto(Vec3 a, Vec3 b, Vec3 c, Vec3 out) {
        // b x c -> tmp
        Vec3 tmp = out; // reuse out as tmp
        float tx = b.y * c.z - b.z * c.y;
        float ty = b.z * c.x - b.x * c.z;
        float tz = b.x * c.y - b.y * c.x;
        // a x tmp -> out
        float x = a.y * tz - a.z * ty;
        float y = a.z * tx - a.x * tz;
        float z = a.x * ty - a.y * tx;
        out.set(x, y, z);
    }

    /** Any perpendicular vector to v (written into out). */
    private static void perpendicularInto(Vec3 v, Vec3 out) {
        if (Math.abs(v.x) < 0.9f) out.set(0, -v.z, v.y);
        else out.set(-v.z, 0, v.x);
    }
}
