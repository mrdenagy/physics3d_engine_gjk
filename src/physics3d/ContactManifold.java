
package physics3d;

import java.util.ArrayList;

/** Persistent contact manifold for a body pair (up to 4 points). */
public final class ContactManifold {
    public final RigidBody A;
    public final RigidBody B;

    public final Vec3 normal = new Vec3(0, 1, 0); // A->B
    public float penetration;

    public float restitution;
    public float friction;

    public final ArrayList<ContactPoint> points = new ArrayList<>(4);
    public int framesSinceUpdate = 0;

    // Reusable temp list to avoid per-update allocations
    private final ArrayList<ContactPoint> tmpNew = new ArrayList<>(4);

    // Scratch vectors to avoid allocations in matching/update/prune
    private final Vec3 tmpWA = new Vec3();
    private final Vec3 tmpWB = new Vec3();
    private final Vec3 tmpMid = new Vec3();
    private final Vec3 tmpLocalA = new Vec3();
    private final Vec3 tmpLocalB = new Vec3();
    private final Vec3 tmpDiff = new Vec3();

    public ContactManifold(RigidBody a, RigidBody b) {
        this.A = a;
        this.B = b;
    }

    public void mixMaterial() {
        restitution = Math.min(A.restitution, B.restitution);
        friction = (float) Math.sqrt(A.friction * B.friction);
    }

    /**
     * Updates this persistent manifold from the narrowphase manifold result.
     *
     * Optimized:
     * - No new ArrayList each call.
     * - No temporary Vec3 allocations from add/sub/mul chains.
     * - Uses allocation-free local/world conversion via RigidBody *Into methods. [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/ContactManifold.java)[1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/ContactManifold.java)
     */
    public void updateFrom(Collision.ManifoldResult r) {
        framesSinceUpdate = 0;
        normal.set(r.normal);
        penetration = r.depth;
        mixMaterial();

        tmpNew.clear();

        int count = Math.min(4, r.pointsA.size());
        // Track which old points were already matched (bitmask over points list)
        int usedMask = 0;

        for (int i = 0; i < count; i++) {
            Vec3 pA = r.pointsA.get(i);
            Vec3 pB = r.pointsB.get(i);

            // cWorld = (pA + pB) * 0.5f (no allocations)
            tmpMid.set(
                    (pA.x + pB.x) * 0.5f,
                    (pA.y + pB.y) * 0.5f,
                    (pA.z + pB.z) * 0.5f
            );

            // local anchors (no allocations)
            A.localPointFromWorldInto(tmpLocalA, pA);
            B.localPointFromWorldInto(tmpLocalB, pB);

            // Find best match among existing points that hasn't been used yet
            int matchIndex = findBestMatchIndex(tmpMid, usedMask);
            if (matchIndex >= 0) {
                usedMask |= (1 << matchIndex);
                ContactPoint match = points.get(matchIndex);

                match.localA.set(tmpLocalA);
                match.localB.set(tmpLocalB);
                match.worldPoint.set(tmpMid);

                tmpNew.add(match);
            } else {
                // No good match: create a new point (rare once things settle)
                ContactPoint cp = new ContactPoint();
                cp.localA.set(tmpLocalA);
                cp.localB.set(tmpLocalB);
                cp.worldPoint.set(tmpMid);
                cp.resetImpulses();
                tmpNew.add(cp);
            }
        }

        points.clear();
        points.addAll(tmpNew);

        pruneStalePoints();
    }

    /**
     * Returns the index of the best matching existing point (or -1 if none).
     * Uses midpoint distance with a fixed threshold. No allocations. [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/ContactManifold.java)
     */
    private int findBestMatchIndex(Vec3 worldPoint, int usedMask) {
        int matchIndex = -1;
        float bestDist2 = Float.POSITIVE_INFINITY;

        for (int i = 0; i < points.size(); i++) {
            if ((usedMask & (1 << i)) != 0) continue;

            ContactPoint cp = points.get(i);

            // wa = A.worldPointFromLocal(cp.localA), wb = B.worldPointFromLocal(cp.localB)
            A.worldPointFromLocalInto(tmpWA, cp.localA);
            B.worldPointFromLocalInto(tmpWB, cp.localB);

            // w = (wa + wb) * 0.5 (no alloc)
            tmpMid.set(
                    (tmpWA.x + tmpWB.x) * 0.5f,
                    (tmpWA.y + tmpWB.y) * 0.5f,
                    (tmpWA.z + tmpWB.z) * 0.5f
            );

            // diff = w - worldPoint
            tmpDiff.set(
                    tmpMid.x - worldPoint.x,
                    tmpMid.y - worldPoint.y,
                    tmpMid.z - worldPoint.z
            );

            float d2 = tmpDiff.len2();
            if (d2 < bestDist2) {
                bestDist2 = d2;
                matchIndex = i;
            }
        }

        float threshold = 0.06f;
        return (matchIndex >= 0 && bestDist2 < threshold * threshold) ? matchIndex : -1;
    }

    /**
     * Recomputes world points for stored contacts and removes stale ones.
     * No allocations. [1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/ContactManifold.java)[1](https://stluciepublicschools-my.sharepoint.com/personal/mde0701_stlucieschools_org/Documents/Microsoft%20Copilot%20Chat%20Files/ContactManifold.java)
     */
    private void pruneStalePoints() {
        for (int i = points.size() - 1; i >= 0; i--) {
            ContactPoint cp = points.get(i);

            A.worldPointFromLocalInto(tmpWA, cp.localA);
            B.worldPointFromLocalInto(tmpWB, cp.localB);

            // cp.worldPoint = midpoint
            cp.worldPoint.set(
                    (tmpWA.x + tmpWB.x) * 0.5f,
                    (tmpWA.y + tmpWB.y) * 0.5f,
                    (tmpWA.z + tmpWB.z) * 0.5f
            );

            // sep = |wa - wb|
            tmpDiff.set(
                    tmpWA.x - tmpWB.x,
                    tmpWA.y - tmpWB.y,
                    tmpWA.z - tmpWB.z
            );
            float sep = (float) Math.sqrt(tmpDiff.len2());

            if (sep > 0.25f) {
                points.remove(i);
            }
        }
    }
}
