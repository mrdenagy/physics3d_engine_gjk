
package physics3d;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;

/**
 * 3D spatial hash broadphase with:
 * - incremental body updates (cell occupancy bounds tracking)
 * - AABB query without per-query allocations
 * - cached pair list (recomputed only when cell occupancy changes)
 *
 * Near-zero GC:
 * - No int[] bounds allocations (uses reusable Bounds objects).
 * - No HashSet<Long> allocations in computePairs (uses reusable LongHashSet).
 * - Pair objects are reused from a pool (grows only if a new maximum is reached).
 * - Avoids computeIfAbsent lambda in tight loops.
 *
 * AABB-cached correctness:
 * - Calls b.updateAABB() before reading b.aabb() in insert/update.
 */
public final class SpatialHashBroadphase3D {
    private final float cellSize;

    private final HashMap<Long, ArrayList<RigidBody>> grid = new HashMap<>(2048);

    /** Store body cell bounds: reusable Bounds object per body. */
    private final IdentityHashMap<RigidBody, Bounds> bodyBounds = new IdentityHashMap<>();

    /** Pair cache */
    private final ArrayList<Pair> cachedPairs = new ArrayList<>(256);
    private boolean pairsDirty = true;

    /** Reusable "seen" set for deduping pairs without allocating HashSet<Long>. */
    private final LongHashSet seenPairs = new LongHashSet(1024);

    /** Reusable Pair pool (grow-only). */
    private final ArrayList<Pair> pairPool = new ArrayList<>(256);

    /** Reusable query bounds (avoid allocations per query). */
    private final Bounds queryBounds = new Bounds();

    public SpatialHashBroadphase3D(float cellSize) {
        this.cellSize = cellSize;
    }

    public void clear() {
        grid.clear();
        bodyBounds.clear();
        cachedPairs.clear();
        pairPool.clear();
        seenPairs.clear();
        pairsDirty = true;
    }

    /** Inserts all bodies from scratch. */
    public void rebuild(List<RigidBody> bodies) {
        clear();
        for (RigidBody b : bodies) insert(b);
        pairsDirty = true;
    }

    public void insert(RigidBody b) {
        // Ensure cached AABB is current (allocation-free).
        b.updateAABB();

        Bounds bounds = new Bounds();
        bounds.computeFromAABB(b.aabb(), cellSize);
        bodyBounds.put(b, bounds);

        addToCells(b, bounds);
        pairsDirty = true;
    }

    /** Incrementally update a body's cell occupancy after it moves. */
    public void update(RigidBody b) {
        Bounds old = bodyBounds.get(b);
        if (old == null) {
            insert(b);
            return;
        }

        // Ensure cached AABB is current (allocation-free).
        b.updateAABB();

        // Compute new bounds into reusable temp and compare.
        queryBounds.computeFromAABB(b.aabb(), cellSize);
        if (old.sameAs(queryBounds)) return;

        removeFromCells(b, old);
        old.set(queryBounds);
        addToCells(b, old);
        pairsDirty = true;
    }

    /** Query bodies overlapping an AABB (fills outUnique, no allocations). */
    public void queryAABBFill(AABB aabb,
                              ArrayList<RigidBody> outUnique,
                              IdentityHashMap<RigidBody, Boolean> seen) {
        outUnique.clear();
        seen.clear();

        queryBounds.computeFromAABB(aabb, cellSize);

        for (int z = queryBounds.z0; z <= queryBounds.z1; z++) {
            for (int y = queryBounds.y0; y <= queryBounds.y1; y++) {
                for (int x = queryBounds.x0; x <= queryBounds.x1; x++) {
                    ArrayList<RigidBody> cell = grid.get(pack(x, y, z));
                    if (cell == null) continue;
                    for (int i = 0; i < cell.size(); i++) {
                        RigidBody rb = cell.get(i);
                        if (seen.put(rb, Boolean.TRUE) == null) outUnique.add(rb);
                    }
                }
            }
        }
    }

    /** Convenience: query bodies overlapping an AABB (allocates) — avoid in near-zero GC code. */
    public ArrayList<RigidBody> queryAABB(AABB aabb) {
        ArrayList<RigidBody> out = new ArrayList<>();
        IdentityHashMap<RigidBody, Boolean> seen = new IdentityHashMap<>();
        queryAABBFill(aabb, out, seen);
        return out;
    }

    /**
     * Compute potentially colliding pairs in the current grid.
     * Pair cache optimization: only recompute when pairsDirty.
     *
     * Near-zero GC: reuses Pair objects and a reusable long set for dedupe.
     */
    public ArrayList<Pair> computePairs() {
        if (!pairsDirty) return cachedPairs;

        cachedPairs.clear();
        seenPairs.clear();

        int outIndex = 0;

        for (ArrayList<RigidBody> cell : grid.values()) {
            int n = cell.size();
            for (int i = 0; i < n; i++) {
                RigidBody a = cell.get(i);
                for (int j = i + 1; j < n; j++) {
                    RigidBody b = cell.get(j);
                    if (a.isStatic() && b.isStatic()) continue;

                    long h = pairHash(a, b);
                    if (!seenPairs.add(h)) continue;

                    Pair p;
                    if (outIndex < pairPool.size()) {
                        p = pairPool.get(outIndex);
                    } else {
                        p = new Pair();
                        pairPool.add(p);
                    }
                    p.a = a;
                    p.b = b;

                    if (outIndex < cachedPairs.size()) cachedPairs.set(outIndex, p);
                    else cachedPairs.add(p);

                    outIndex++;
                }
            }
        }

        // Trim logical size without freeing capacity
        while (cachedPairs.size() > outIndex) cachedPairs.remove(cachedPairs.size() - 1);

        pairsDirty = false;
        return cachedPairs;
    }

    public static final class Pair {
        public RigidBody a, b;
        public Pair() {}
        public Pair(RigidBody a, RigidBody b) { this.a = a; this.b = b; }
    }

    public static long pairHash(RigidBody a, RigidBody b) {
        int ha = System.identityHashCode(a);
        int hb = System.identityHashCode(b);
        int lo = Math.min(ha, hb), hi = Math.max(ha, hb);
        return (((long) lo) << 32) ^ (hi & 0xffffffffL);
    }

    // --- internal helpers ---

    private long pack(int x, int y, int z) {
        long k = 1469598103934665603L;
        k ^= x; k *= 1099511628211L;
        k ^= y; k *= 1099511628211L;
        k ^= z; k *= 1099511628211L;
        return k;
    }

    private void addToCells(RigidBody b, Bounds bounds) {
        for (int z = bounds.z0; z <= bounds.z1; z++) {
            for (int y = bounds.y0; y <= bounds.y1; y++) {
                for (int x = bounds.x0; x <= bounds.x1; x++) {
                    long k = pack(x, y, z);
                    ArrayList<RigidBody> list = grid.get(k);
                    if (list == null) {
                        list = new ArrayList<>(4);
                        grid.put(k, list);
                    }
                    list.add(b);
                }
            }
        }
    }

    private void removeFromCells(RigidBody b, Bounds bounds) {
        for (int z = bounds.z0; z <= bounds.z1; z++) {
            for (int y = bounds.y0; y <= bounds.y1; y++) {
                for (int x = bounds.x0; x <= bounds.x1; x++) {
                    long k = pack(x, y, z);
                    ArrayList<RigidBody> list = grid.get(k);
                    if (list == null) continue;

                    for (int i = list.size() - 1; i >= 0; i--) {
                        if (list.get(i) == b) { list.remove(i); break; }
                    }
                    if (list.isEmpty()) grid.remove(k);
                }
            }
        }
    }

    /** Reusable integer bounds (no int[] allocations). */
    private static final class Bounds {
        int x0, y0, z0, x1, y1, z1;

        void set(Bounds o) {
            x0 = o.x0; y0 = o.y0; z0 = o.z0;
            x1 = o.x1; y1 = o.y1; z1 = o.z1;
        }

        boolean sameAs(Bounds o) {
            return x0 == o.x0 && y0 == o.y0 && z0 == o.z0 &&
                   x1 == o.x1 && y1 == o.y1 && z1 == o.z1;
        }

        void computeFromAABB(AABB a, float cellSize) {
            x0 = (int) Math.floor(a.min.x / cellSize);
            y0 = (int) Math.floor(a.min.y / cellSize);
            z0 = (int) Math.floor(a.min.z / cellSize);

            x1 = (int) Math.floor(a.max.x / cellSize);
            y1 = (int) Math.floor(a.max.y / cellSize);
            z1 = (int) Math.floor(a.max.z / cellSize);
        }
    }

    /**
     * Reusable primitive long hash set (open addressing).
     * - No boxing
     * - No per-frame allocation
     * - Grows only when needed (rare after warm-up)
     */
    private static final class LongHashSet {
        private long[] keys;
        private boolean[] used;
        private int mask;
        private int size;
        private int threshold;

        LongHashSet(int initialCapacity) {
            int cap = 1;
            while (cap < initialCapacity) cap <<= 1;
            keys = new long[cap];
            used = new boolean[cap];
            mask = cap - 1;
            threshold = (int) (cap * 0.7f);
        }

        void clear() {
            if (size == 0) return;
            for (int i = 0; i < used.length; i++) used[i] = false;
            size = 0;
        }

        boolean add(long k) {
            if (size >= threshold) rehash();

            int idx = mix64(k) & mask;
            while (used[idx]) {
                if (keys[idx] == k) return false;
                idx = (idx + 1) & mask;
            }
            used[idx] = true;
            keys[idx] = k;
            size++;
            return true;
        }

        private void rehash() {
            long[] oldK = keys;
            boolean[] oldU = used;

            int newCap = oldK.length << 1;
            keys = new long[newCap];
            used = new boolean[newCap];
            mask = newCap - 1;
            threshold = (int) (newCap * 0.7f);
            size = 0;

            for (int i = 0; i < oldK.length; i++) {
                if (!oldU[i]) continue;
                add(oldK[i]);
            }
        }

        private static int mix64(long z) {
            z ^= (z >>> 33);
            z *= 0xff51afd7ed558ccdL;
            z ^= (z >>> 33);
            z *= 0xc4ceb9fe1a85ec53L;
            z ^= (z >>> 33);
            return (int) z;
        }
    }
}
