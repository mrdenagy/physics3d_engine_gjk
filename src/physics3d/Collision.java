
package physics3d;

import java.util.ArrayList;

public final class Collision {
    private Collision() {}

    public static final class ManifoldResult {
        public boolean hit;
        public final Vec3 normal = new Vec3(0, 1, 0); // A->B
        public float depth;

        // Exposed as lists for compatibility with existing code (ContactManifold.updateFrom, etc.)
        public final ArrayList<Vec3> pointsA = new ArrayList<>(4);
        public final ArrayList<Vec3> pointsB = new ArrayList<>(4);

        // Internal reusable pools so we can clear() without re-allocating Vec3s.
        private final Vec3[] poolA = new Vec3[] { new Vec3(), new Vec3(), new Vec3(), new Vec3() };
        private final Vec3[] poolB = new Vec3[] { new Vec3(), new Vec3(), new Vec3(), new Vec3() };

        public ManifoldResult() {}

        public void clear() {
            hit = false;
            depth = 0f;
            normal.set(0, 1, 0);
            pointsA.clear();
            pointsB.clear();
        }

        /** Ensures lists contain exactly count points, reusing the same Vec3 instances each time. */
        public void setPointCount(int count) {
            pointsA.clear();
            pointsB.clear();
            count = Math.max(0, Math.min(4, count));
            for (int i = 0; i < count; i++) {
                pointsA.add(poolA[i]);
                pointsB.add(poolB[i]);
            }
        }
    }

    /** Backwards-compatible allocating version (kept). [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/World3D.java) */
    public static ManifoldResult detectManifold(RigidBody A, RigidBody B) {
        ManifoldResult r = new ManifoldResult();
        detectManifoldInto(A, B, r);
        return r;
    }

    /**
     * Allocation-friendly version: fills the provided output manifold without creating
     * a new ManifoldResult each call.
     */
    public static void detectManifoldInto(RigidBody A, RigidBody B, ManifoldResult out) {
        out.clear();

        // Box-box special case (SAT + clipping gives up to 4 contacts)
        if (A.shape.type() == Shape.Type.BOX && B.shape.type() == Shape.Type.BOX) {
            BoxBoxManifold.ManifoldResult bb =
                    BoxBoxManifold.collide(A, B, (BoxShape) A.shape, (BoxShape) B.shape);

            out.hit = bb.hit;
            if (!out.hit) return;

            out.normal.set(bb.normal);
            out.depth = bb.penetration;

            int count = Math.min(4, bb.pointsA.size());
            out.setPointCount(count);
            for (int i = 0; i < count; i++) {
                out.pointsA.get(i).set(bb.pointsA.get(i));
                out.pointsB.get(i).set(bb.pointsB.get(i));
            }
            return;
        }

        // General convex path: GJK -> EPA
        GJK.Result gjk = GJK.intersect(A, B);
        if (!gjk.intersect) {
            out.hit = false;
            return;
        }

        EPA.ContactInfo info = EPA.expand(A, B, gjk.simplex);

        out.hit = true;
        out.normal.set(info.normal);
        out.depth = info.depth;

        out.setPointCount(1);
        out.pointsA.get(0).set(info.pointA);
        out.pointsB.get(0).set(info.pointB);
    }
}
