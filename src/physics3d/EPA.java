package physics3d;

import java.util.ArrayList;
import java.util.List;

/**
 * EPA (Expanding Polytope Algorithm) for penetration depth/normal after GJK reports intersection.
 *
 * Optimized:
 * - Uses Minkowski.supportInto(...) (allocation-free support mapping) instead of Minkowski.support(...). [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/EPA.java)
 * - Avoids allocating Vec3 temporaries in hot loops by using explicit float math + Vec3.set.
 * - Pools Faces and SupportPoints using ThreadLocal Workspace.
 * - Avoids allocating horizon edges as int[].
 */
public final class EPA {
    private EPA() {}

    public static final class ContactInfo {
        public Vec3 normal = new Vec3(0, 1, 0);   // from A -> B
        public float depth = 0f;
        public Vec3 pointA = new Vec3();
        public Vec3 pointB = new Vec3();
        public Vec3 contactPoint = new Vec3();    // midpoint
    }

    private static final class Face {
        int a, b, c;     // indices into pts
        final Vec3 n = new Vec3(); // outward normal
        float d;         // distance to origin along normal

        void set(int a, int b, int c) { this.a = a; this.b = b; this.c = c; }
    }

    private static final class Workspace {
        // Pools / reusable lists
        final ArrayList<Minkowski.SupportPoint> pts = new ArrayList<>(64);
        final ArrayList<Face> faces = new ArrayList<>(128);

        // Pools
        final Minkowski.SupportPoint[] spPool = new Minkowski.SupportPoint[128];
        int spUsed = 0;

        final Face[] facePool = new Face[256];
        int faceUsed = 0;

        // Horizon edge buffer (pairs of ints); cancellation removes reverse edges
        final int[] edgeA = new int[1024];
        final int[] edgeB = new int[1024];
        int edgeCount = 0;

        // Scratch vectors
        final Vec3 dir = new Vec3(1, 0, 0);
        final Vec3 tmp = new Vec3();
        final Vec3 tmp2 = new Vec3();
        final Vec3 tmp3 = new Vec3();
        final Vec3 p = new Vec3();     // closest point on face plane to origin
        final Vec3 v0 = new Vec3();
        final Vec3 v1 = new Vec3();
        final Vec3 v2 = new Vec3();

        Workspace() {
            for (int i = 0; i < spPool.length; i++) spPool[i] = new Minkowski.SupportPoint();
            for (int i = 0; i < facePool.length; i++) facePool[i] = new Face();
        }

        void reset() {
            pts.clear();
            faces.clear();
            spUsed = 0;
            faceUsed = 0;
            edgeCount = 0;
        }

        Minkowski.SupportPoint newSupport() {
            if (spUsed >= spPool.length) {
                // Fallback (rare): grow pool by direct allocation if scene is extreme
                return new Minkowski.SupportPoint();
            }
            return spPool[spUsed++];
        }

        Face newFace() {
            if (faceUsed >= facePool.length) {
                // Fallback (rare)
                return new Face();
            }
            return facePool[faceUsed++];
        }
    }

    private static final ThreadLocal<Workspace> TLS = ThreadLocal.withInitial(Workspace::new);

    public static ContactInfo expand(RigidBody A, RigidBody B, List<Minkowski.SupportPoint> simplex) {
        Workspace ws = TLS.get();
        ws.reset();

        // Copy simplex into EPA-owned support points so we don't mutate GJK's simplex objects.
        for (int i = 0; i < simplex.size(); i++) {
            Minkowski.SupportPoint src = simplex.get(i);
            Minkowski.SupportPoint dst = ws.newSupport();
            copySupport(dst, src);
            ws.pts.add(dst);
        }

        // Ensure tetrahedron (4 points)
        while (ws.pts.size() < 4) {
            Minkowski.SupportPoint base = ws.pts.get(0);
            // dir = normalize(base.v) with robustness
            normalizeInto(base.v, ws.dir);

            if (ws.dir.len2() < 1e-8f) ws.dir.set(1, 0, 0);

            Minkowski.SupportPoint p = ws.newSupport();
            Minkowski.supportInto(A, B, ws.dir, p);
            ws.pts.add(p);
        }

        // Initial tetra faces
        addFace(ws, 0, 1, 2);
        addFace(ws, 0, 3, 1);
        addFace(ws, 0, 2, 3);
        addFace(ws, 1, 3, 2);

        // Expand
        for (int iter = 0; iter < 48; iter++) {
            Face f = closestFace(ws.faces);
            // dir is face normal
            ws.dir.set(f.n);

            Minkowski.SupportPoint newPt = ws.newSupport();
            Minkowski.supportInto(A, B, ws.dir, newPt);

            float dist = dot(newPt.v, ws.dir);
            if (dist - f.d < 1e-4f) {
                return buildContactInfo(ws, f);
            }

            int newIndex = ws.pts.size();
            ws.pts.add(newPt);

            // Build horizon edges by removing faces visible from new point
            ws.edgeCount = 0;

            for (int i = ws.faces.size() - 1; i >= 0; i--) {
                Face face = ws.faces.get(i);

                Vec3 va = ws.pts.get(face.a).v;
                // visible if face.n · (p - va) > 0
                subInto(newPt.v, va, ws.tmp);
                if (dot(face.n, ws.tmp) > 0f) {
                    addEdge(ws, face.a, face.b);
                    addEdge(ws, face.b, face.c);
                    addEdge(ws, face.c, face.a);
                    ws.faces.remove(i);
                }
            }

            // Rebuild polytope with new faces from horizon edges
            for (int e = 0; e < ws.edgeCount; e++) {
                addFace(ws, ws.edgeA[e], ws.edgeB[e], newIndex);
            }
        }

        // Fallback: return best face
        return buildContactInfo(ws, closestFace(ws.faces));
    }

    // ----------------- Face / Polytope helpers -----------------

    private static void addFace(Workspace ws, int a, int b, int c) {
        Face f = ws.newFace();
        f.set(a, b, c);
        computeFacePlane(ws, f);

        // Ensure outward normal so d is positive
        if (f.d < 0f) {
            int tmp = f.b; f.b = f.c; f.c = tmp;
            computeFacePlane(ws, f);
        }

        ws.faces.add(f);
    }

    private static void computeFacePlane(Workspace ws, Face f) {
        Vec3 A = ws.pts.get(f.a).v;
        Vec3 B = ws.pts.get(f.b).v;
        Vec3 C = ws.pts.get(f.c).v;

        // n = normalize((B-A) x (C-A))
        subInto(B, A, ws.tmp);
        subInto(C, A, ws.tmp2);
        crossInto(ws.tmp, ws.tmp2, f.n);
        normalizeSelf(f.n);

        f.d = dot(f.n, A);
    }

    private static Face closestFace(ArrayList<Face> faces) {
        Face best = faces.get(0);
        for (int i = 1; i < faces.size(); i++) {
            Face f = faces.get(i);
            if (f.d < best.d) best = f;
        }
        return best;
    }

    /**
     * Add edge (a->b). If reverse edge already exists (b->a), remove it (internal edge).
     * Uses Workspace edge buffers to avoid allocating int[].
     */
    private static void addEdge(Workspace ws, int a, int b) {
        // cancel reverse edge if present
        for (int i = 0; i < ws.edgeCount; i++) {
            if (ws.edgeA[i] == b && ws.edgeB[i] == a) {
                // remove by swap-with-last
                int last = ws.edgeCount - 1;
                ws.edgeA[i] = ws.edgeA[last];
                ws.edgeB[i] = ws.edgeB[last];
                ws.edgeCount--;
                return;
            }
        }
        if (ws.edgeCount >= ws.edgeA.length) return; // safety
        ws.edgeA[ws.edgeCount] = a;
        ws.edgeB[ws.edgeCount] = b;
        ws.edgeCount++;
    }

    // ----------------- Contact construction -----------------

    private static ContactInfo buildContactInfo(Workspace ws, Face f) {
        Vec3 A = ws.pts.get(f.a).v;
        Vec3 B = ws.pts.get(f.b).v;
        Vec3 C = ws.pts.get(f.c).v;

        Vec3 n = f.n;

        // p = n * d (closest point to origin on plane along normal)
        ws.p.set(n.x * f.d, n.y * f.d, n.z * f.d);

        // Compute barycentric of p in triangle ABC (in Minkowski space)
        // v0 = B-A, v1 = C-A, v2 = p-A
        subInto(B, A, ws.v0);
        subInto(C, A, ws.v1);
        subInto(ws.p, A, ws.v2);

        float d00 = dot(ws.v0, ws.v0);
        float d01 = dot(ws.v0, ws.v1);
        float d11 = dot(ws.v1, ws.v1);
        float d20 = dot(ws.v2, ws.v0);
        float d21 = dot(ws.v2, ws.v1);

        float denom = d00 * d11 - d01 * d01;
        float invDen = (Math.abs(denom) < 1e-8f) ? 1f : (1f / denom);

        float v = (d11 * d20 - d01 * d21) * invDen;
        float w = (d00 * d21 - d01 * d20) * invDen;
        float u = 1f - v - w;

        Minkowski.SupportPoint spA = ws.pts.get(f.a);
        Minkowski.SupportPoint spB = ws.pts.get(f.b);
        Minkowski.SupportPoint spC = ws.pts.get(f.c);

        // pointOnA = u*aA + v*aB + w*aC
        Vec3 pointOnA = ws.tmp;   // reuse scratch as output holder
        Vec3 pointOnB = ws.tmp2;  // reuse scratch as output holder

        pointOnA.set(
                spA.a.x * u + spB.a.x * v + spC.a.x * w,
                spA.a.y * u + spB.a.y * v + spC.a.y * w,
                spA.a.z * u + spB.a.z * v + spC.a.z * w
        );
        pointOnB.set(
                spA.b.x * u + spB.b.x * v + spC.b.x * w,
                spA.b.y * u + spB.b.y * v + spC.b.y * w,
                spA.b.z * u + spB.b.z * v + spC.b.z * w
        );

        ContactInfo ci = new ContactInfo();
        ci.normal.set(n);
        ci.depth = f.d;
        ci.pointA.set(pointOnA);
        ci.pointB.set(pointOnB);

        // contactPoint = (pointA + pointB) * 0.5
        ci.contactPoint.set(
                (pointOnA.x + pointOnB.x) * 0.5f,
                (pointOnA.y + pointOnB.y) * 0.5f,
                (pointOnA.z + pointOnB.z) * 0.5f
        );
        return ci;
    }

    // ----------------- Small vector helpers (avoid Vec3 allocations) -----------------

    private static void copySupport(Minkowski.SupportPoint dst, Minkowski.SupportPoint src) {
        dst.v.set(src.v);
        dst.a.set(src.a);
        dst.b.set(src.b);
    }

    private static float dot(Vec3 a, Vec3 b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    private static void subInto(Vec3 a, Vec3 b, Vec3 out) {
        out.set(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    private static void crossInto(Vec3 a, Vec3 b, Vec3 out) {
        out.set(
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x
        );
    }

    private static void normalizeInto(Vec3 v, Vec3 out) {
        float x = v.x, y = v.y, z = v.z;
        float len2 = x * x + y * y + z * z;
        if (len2 < 1e-12f) {
            out.set(0, 0, 0);
            return;
        }
        float inv = 1f / (float) Math.sqrt(len2);
        out.set(x * inv, y * inv, z * inv);
    }

    private static void normalizeSelf(Vec3 v) {
        float len2 = v.x * v.x + v.y * v.y + v.z * v.z;
        if (len2 < 1e-12f) return;
        float inv = 1f / (float) Math.sqrt(len2);
        v.set(v.x * inv, v.y * inv, v.z * inv);
    }
}
