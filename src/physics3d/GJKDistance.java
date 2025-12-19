package physics3d;

/**
 * Allocation-reduced GJK distance algorithm for convex shapes using support mappings.
 *
 * This version calls Shape.supportInto(...) to avoid Vec3 allocations from support mapping.
 */
public final class GJKDistance {
    private GJKDistance() {}

    public static final class Result {
        public boolean intersect;
        public float distance;
        public final Vec3 normal = new Vec3(1,0,0); // from A to B
        public final Vec3 pointA = new Vec3();
        public final Vec3 pointB = new Vec3();
    }

    private static final class Vertex {
        final Vec3 v = new Vec3();
        final Vec3 a = new Vec3();
        final Vec3 b = new Vec3();
    }

    private static final class Workspace {
        final Vertex[] s = new Vertex[]{ new Vertex(), new Vertex(), new Vertex(), new Vertex() };
        int n = 0;

        final Vec3 dir = new Vec3(1,0,0);

        // scratch vectors to keep support mapping allocation-free
        final Vec3 negDir = new Vec3();
        final Vec3 supA = new Vec3();
        final Vec3 supB = new Vec3();

        // barycentric
        float u,v,w;
    }

    private static final ThreadLocal<Workspace> TLS = ThreadLocal.withInitial(Workspace::new);

    public static Result distance(RigidBody A, RigidBody B){
        Workspace ws = TLS.get();
        ws.n = 0;

        Result res = new Result();

        ws.dir.set(B.position.x - A.position.x, B.position.y - A.position.y, B.position.z - A.position.z);
        if (ws.dir.len2() < 1e-10f) ws.dir.set(1,0,0);

        supportInto(A, B, ws, ws.dir, ws.s[0]);
        ws.n = 1;

        ws.dir.set(-ws.s[0].v.x, -ws.s[0].v.y, -ws.s[0].v.z);

        float prevDist2 = Float.POSITIVE_INFINITY;

        for (int iter=0; iter<32; iter++){
            float dir2 = ws.dir.len2();
            if (dir2 < 1e-14f){
                res.intersect = true;
                res.distance = 0f;
                res.normal.set(1,0,0);
                return res;
            }

            Vertex vNew = ws.s[Math.min(ws.n, 3)];
            supportInto(A, B, ws, ws.dir, vNew);

            float invLen = 1f / (float)Math.sqrt(dir2);
            float dirHx = ws.dir.x * invLen;
            float dirHy = ws.dir.y * invLen;
            float dirHz = ws.dir.z * invLen;
            float proj = vNew.v.x*dirHx + vNew.v.y*dirHy + vNew.v.z*dirHz;
            if (proj > 0f && proj*proj >= prevDist2 - 1e-8f) {
                break;
            }

            if (ws.n == 1){
                swap(ws.s[1], ws.s[0]);
                copyVertex(ws.s[0], vNew);
                ws.n = 2;
            } else if (ws.n == 2){
                swap(ws.s[2], ws.s[1]);
                swap(ws.s[1], ws.s[0]);
                copyVertex(ws.s[0], vNew);
                ws.n = 3;
            } else if (ws.n == 3){
                swap(ws.s[3], ws.s[2]);
                swap(ws.s[2], ws.s[1]);
                swap(ws.s[1], ws.s[0]);
                copyVertex(ws.s[0], vNew);
                ws.n = 4;
            } else {
                swap(ws.s[3], ws.s[2]);
                swap(ws.s[2], ws.s[1]);
                swap(ws.s[1], ws.s[0]);
                copyVertex(ws.s[0], vNew);
            }

            boolean contains = closestToOrigin(ws, res);
            if (contains){
                res.intersect = true;
                res.distance = 0f;
                return res;
            }

            prevDist2 = ws.dir.len2();
        }

        closestToOrigin(ws, res);
        res.intersect = false;
        float d2 = ws.dir.len2();
        res.distance = (float)Math.sqrt(Math.max(0f, d2));
        if (res.distance > 1e-8f){
            float inv = 1f/res.distance;
            res.normal.set(ws.dir.x*inv, ws.dir.y*inv, ws.dir.z*inv);
        }
        return res;
    }

    private static void supportInto(RigidBody A, RigidBody B, Workspace ws, Vec3 dir, Vertex out){
        A.shape.supportInto(ws.supA, dir, A.position, A.orientation);
        ws.negDir.set(-dir.x, -dir.y, -dir.z);
        B.shape.supportInto(ws.supB, ws.negDir, B.position, B.orientation);

        out.a.set(ws.supA);
        out.b.set(ws.supB);
        out.v.set(ws.supA.x - ws.supB.x, ws.supA.y - ws.supB.y, ws.supA.z - ws.supB.z);
    }

    private static void copyVertex(Vertex dst, Vertex src){
        dst.v.set(src.v);
        dst.a.set(src.a);
        dst.b.set(src.b);
    }

    private static void swap(Vertex a, Vertex b){
        float tx=a.v.x, ty=a.v.y, tz=a.v.z;
        a.v.set(b.v);
        b.v.set(tx,ty,tz);

        tx=a.a.x; ty=a.a.y; tz=a.a.z;
        a.a.set(b.a);
        b.a.set(tx,ty,tz);

        tx=a.b.x; ty=a.b.y; tz=a.b.z;
        a.b.set(b.b);
        b.b.set(tx,ty,tz);
    }

    private static boolean closestToOrigin(Workspace ws, Result res){
        if (ws.n == 1){
            Vertex A = ws.s[0];
            ws.dir.set(-A.v.x, -A.v.y, -A.v.z);
            res.pointA.set(A.a);
            res.pointB.set(A.b);
            res.normal.set(A.v.normalized());
            return ws.dir.len2() < 1e-12f;
        }

        if (ws.n == 2){
            Vertex A = ws.s[0];
            Vertex B = ws.s[1];

            float abx = B.v.x - A.v.x;
            float aby = B.v.y - A.v.y;
            float abz = B.v.z - A.v.z;
            float ab2 = abx*abx + aby*aby + abz*abz;

            float t = 0f;
            if (ab2 > 1e-12f){
                t = -(A.v.x*abx + A.v.y*aby + A.v.z*abz) / ab2;
                t = MathUtil.clamp(t, 0f, 1f);
            }

            float cx = A.v.x + abx*t;
            float cy = A.v.y + aby*t;
            float cz = A.v.z + abz*t;
            ws.dir.set(-cx, -cy, -cz);

            float u = 1f - t;
            float v = t;

            res.pointA.set(
                    A.a.x*u + B.a.x*v,
                    A.a.y*u + B.a.y*v,
                    A.a.z*u + B.a.z*v
            );
            res.pointB.set(
                    A.b.x*u + B.b.x*v,
                    A.b.y*u + B.b.y*v,
                    A.b.z*u + B.b.z*v
            );

            if (t <= 1e-5f){
                ws.n = 1;
            } else if (t >= 1f - 1e-5f){
                copyVertex(ws.s[0], B);
                ws.n = 1;
            }

            res.normal.set(cx, cy, cz);
            float l2 = cx*cx + cy*cy + cz*cz;
            if (l2 > 1e-20f){
                float inv = 1f/(float)Math.sqrt(l2);
                res.normal.mulLocal(inv);
            }
            return l2 < 1e-12f;
        }

        if (ws.n == 3){
            Vertex A = ws.s[0];
            Vertex B = ws.s[1];
            Vertex C = ws.s[2];

            triClosestToOrigin(ws, A.v, B.v, C.v);
            float u = ws.u, v = ws.v, w = ws.w;

            float cx = A.v.x*u + B.v.x*v + C.v.x*w;
            float cy = A.v.y*u + B.v.y*v + C.v.y*w;
            float cz = A.v.z*u + B.v.z*v + C.v.z*w;
            float l2 = cx*cx + cy*cy + cz*cz;
            ws.dir.set(-cx, -cy, -cz);

            res.pointA.set(
                    A.a.x*u + B.a.x*v + C.a.x*w,
                    A.a.y*u + B.a.y*v + C.a.y*w,
                    A.a.z*u + B.a.z*v + C.a.z*w
            );
            res.pointB.set(
                    A.b.x*u + B.b.x*v + C.b.x*w,
                    A.b.y*u + B.b.y*v + C.b.y*w,
                    A.b.z*u + B.b.z*v + C.b.z*w
            );

            final float eps = 1e-5f;
            if (u < eps){
                copyVertex(ws.s[0], B);
                copyVertex(ws.s[1], C);
                ws.n = 2;
            } else if (v < eps){
                copyVertex(ws.s[1], C);
                ws.n = 2;
            } else if (w < eps){
                ws.n = 2;
            }

            res.normal.set(cx, cy, cz);
            if (l2 > 1e-20f){
                float inv = 1f/(float)Math.sqrt(l2);
                res.normal.mulLocal(inv);
            }
            return l2 < 1e-12f;
        }

        res.pointA.set(ws.s[0].a);
        res.pointB.set(ws.s[0].b);
        ws.dir.set(0,0,0);
        res.normal.set(1,0,0);
        return true;
    }

    private static void triClosestToOrigin(Workspace ws, Vec3 A, Vec3 B, Vec3 C){
        float abx = B.x - A.x, aby = B.y - A.y, abz = B.z - A.z;
        float acx = C.x - A.x, acy = C.y - A.y, acz = C.z - A.z;

        float aox = -A.x, aoy = -A.y, aoz = -A.z;

        float d1 = abx*aox + aby*aoy + abz*aoz;
        float d2 = acx*aox + acy*aoy + acz*aoz;
        if (d1 <= 0f && d2 <= 0f){
            ws.u = 1f; ws.v = 0f; ws.w = 0f;
            return;
        }

        float box = -B.x, boy = -B.y, boz = -B.z;
        float d3 = abx*box + aby*boy + abz*boz;
        float d4 = acx*box + acy*boy + acz*boz;
        if (d3 >= 0f && d4 <= d3){
            ws.u = 0f; ws.v = 1f; ws.w = 0f;
            return;
        }

        float vc = d1*d4 - d3*d2;
        if (vc <= 0f && d1 >= 0f && d3 <= 0f){
            float v = d1 / (d1 - d3);
            ws.u = 1f - v; ws.v = v; ws.w = 0f;
            return;
        }

        float cox = -C.x, coy = -C.y, coz = -C.z;
        float d5 = abx*cox + aby*coy + abz*coz;
        float d6 = acx*cox + acy*coy + acz*coz;
        if (d6 >= 0f && d5 <= d6){
            ws.u = 0f; ws.v = 0f; ws.w = 1f;
            return;
        }

        float vb = d5*d2 - d1*d6;
        if (vb <= 0f && d2 >= 0f && d6 <= 0f){
            float w = d2 / (d2 - d6);
            ws.u = 1f - w; ws.v = 0f; ws.w = w;
            return;
        }

        float va = d3*d6 - d5*d4;
        float d43 = d4 - d3;
        float d56 = d5 - d6;
        if (va <= 0f && d43 >= 0f && d56 >= 0f){
            float w = d43 / (d43 + d56);
            ws.u = 0f; ws.v = 1f - w; ws.w = w;
            return;
        }

        float denom = 1f / (va + vb + vc);
        float v = vb * denom;
        float w = vc * denom;
        float u = 1f - v - w;
        ws.u = u; ws.v = v; ws.w = w;
    }
}
