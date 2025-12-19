
package physics3d;

import java.util.*;

public final class World3D {
    private final ArrayList<RigidBody> bodies = new ArrayList<>();
    private final ArrayList<ContactManifold> activeManifolds = new ArrayList<>();
    private final HashMap<Long, ContactManifold> manifoldCache = new HashMap<>();
    private final ArrayList<Joint> joints = new ArrayList<>();
    private final SpatialHashBroadphase3D broadphase;

    // Perf: track moved bodies so we only update broadphase for those
    private final IdentityHashMap<RigidBody, Boolean> movedSet = new IdentityHashMap<>();
    private final ArrayList<RigidBody> movedList = new ArrayList<>();

    // Perf: cache swept-AABB query results per bullet per slice
    private static final class SweepCacheEntry {
        public final AABB aabb = new AABB();
        public float dt;
        public float px, py, pz, vx, vy, vz;
        public final ArrayList<RigidBody> candidates = new ArrayList<>();
    }
    private final IdentityHashMap<RigidBody, SweepCacheEntry> sweepCache = new IdentityHashMap<>();
    private final IdentityHashMap<RigidBody, Boolean> querySeen = new IdentityHashMap<>();

    // ✅ Reuse collision manifolds in buildContacts (eliminates per-pair allocations)
    private final Collision.ManifoldResult tmpManifold = new Collision.ManifoldResult();
    private final Collision.ManifoldResult tmpSpecManifold = new Collision.ManifoldResult();

    // ✅ Scratch vectors used in speculative contacts
    private final Vec3 scratchRelV = new Vec3();
    private final Vec3 scratchN = new Vec3();

    // ✅ Solver scratch vectors (hot loop: warmstart + velocity solve + positional solve)
    private final Vec3 sWA = new Vec3();
    private final Vec3 sWB = new Vec3();
    private final Vec3 sP = new Vec3();     // contact world point
    private final Vec3 sPn = new Vec3();    // normal impulse vector
    private final Vec3 sPt = new Vec3();    // tangent impulse vector
    private final Vec3 sImp = new Vec3();   // generic impulse
    private final Vec3 sNegImp = new Vec3();// -impulse

    private final Vec3 sRA = new Vec3();
    private final Vec3 sRB = new Vec3();
    private final Vec3 sOmegaXRA = new Vec3();
    private final Vec3 sOmegaXRB = new Vec3();
    private final Vec3 sVA = new Vec3();
    private final Vec3 sVB = new Vec3();
    private final Vec3 sRV = new Vec3();

    private final Vec3 sRaXn = new Vec3();
    private final Vec3 sRbXn = new Vec3();
    private final Vec3 sInvI_RaXn = new Vec3();
    private final Vec3 sInvI_RbXn = new Vec3();
    private final Vec3 sTmpCross = new Vec3();

    private final Vec3 sT = new Vec3(); // tangent
    private final Vec3 sRaXt = new Vec3();
    private final Vec3 sRbXt = new Vec3();
    private final Vec3 sInvI_RaXt = new Vec3();
    private final Vec3 sInvI_RbXt = new Vec3();

    private final Vec3 sCorrection = new Vec3();

    // Scratch for inv inertia multiply
    private final Vec3 sBody = new Vec3();
    private final Vec3 sScaled = new Vec3();

    // Scratch used in force integration (angular acceleration)
    private final Vec3 sAngAcc = new Vec3();

    // ✅ NEW: allocation-free quaternion integration scratch (no new Quat per body per step)
    private final Quat qDQ = new Quat();     // delta rotation
    private final Quat qNext = new Quat();   // next orientation
    private static final float QUAT_NORM_EPS = 1e-8f;

    // ✅ NEW: broadphase built flag so we don’t rebuild every frame when incrementalBroadphase=true
    private boolean broadphaseBuilt = false;

    public Vec3 gravity = new Vec3(0, -9.8f, 0);
    public int velocityIterations = 18;
    public int positionIterations = 2;

    /** Extra position correction iterations applied per CCD substep (TOI polishing). */
    public int substepPositionIterations = 1;

    public int maxStaleFrames = 10;
    public float restitutionThreshold = 1.0f;
    public float baumgarte = 0.20f;
    public float biasSlop = 0.01f;

    public float linearSleepTolerance = 0.05f;
    public float angularSleepTolerance = 0.05f;
    public float timeToSleep = 0.5f;

    // CCD
    public boolean enableCCD = true;
    public boolean ccdUseAngularBound = true;
    public int maxTOISteps = 4;
    public float minSubstep = 1e-4f;

    // Speculative contacts
    public boolean enableSpeculativeContacts = true;
    public float speculativeFactor = 1.0f;
    public float speculativeSlop = 0.02f;
    public boolean speculativeRequireClosingSpeed = true;
    public float speculativeClosingSpeedThreshold = 0.1f;

    // Perf: broadphase incremental updates
    public boolean incrementalBroadphase = true;

    public World3D(float broadphaseCellSize) {
        this.broadphase = new SpatialHashBroadphase3D(broadphaseCellSize);
    }

    public void add(RigidBody b) {
        bodies.add(b);
        // If bodies are added after stepping begins, force a rebuild next step.
        broadphaseBuilt = false;
    }

    public void addJoint(Joint j) { joints.add(j); }

    public ArrayList<RigidBody> bodies() { return bodies; }
    public ArrayList<Joint> joints() { return joints; }

    public void step(float dt) {
        for (ContactManifold m : manifoldCache.values()) m.framesSinceUpdate++;

        // integrate forces -> velocities
        for (RigidBody b : bodies) {
            if (b.isStatic() || !b.awake) continue;

            b.linearVelocity.x += (gravity.x + b.force.x * b.invMass) * dt;
            b.linearVelocity.y += (gravity.y + b.force.y * b.invMass) * dt;
            b.linearVelocity.z += (gravity.z + b.force.z * b.invMass) * dt;

            // angularVelocity += invI_world * torque * dt (allocation-free)
            invInertiaWorldMulInto(sAngAcc, b, b.torque);
            sAngAcc.mulLocal(dt);
            b.angularVelocity.addLocal(sAngAcc);

            // Damping
            float ld = b.linearDamping;
            float ad = b.angularDamping;
            if (ld > 0f) {
                float k = 1f / (1f + ld * dt);
                b.linearVelocity.mulLocal(k);
            }
            if (ad > 0f) {
                float k = 1f / (1f + ad * dt);
                b.angularVelocity.mulLocal(k);
            }

            b.force.set(0, 0, 0);
            b.torque.set(0, 0, 0);
        }

        manifoldCache.entrySet().removeIf(e -> e.getValue().framesSinceUpdate > maxStaleFrames);

        // ✅ Broadphase: rebuild only once (or whenever bodies list changed), then rely on incremental updates
        if (incrementalBroadphase) {
            if (!broadphaseBuilt) {
                broadphase.rebuild(bodies);
                broadphaseBuilt = true;
            }
        } else {
            // Non-incremental mode: rebuild every frame (debug / comparison mode)
            broadphase.rebuild(bodies);
        }

        movedSet.clear();
        movedList.clear();

        float remaining = dt;
        int toiCount = 0;
        boolean hadAnyTOI = false;

        while (remaining > 0f) {
            // Build contacts at current time
            buildContacts(remaining);

            // Warm start
            for (ContactManifold m : activeManifolds) warmStartContacts(m);
            for (Joint j : joints) j.warmStart();

            // Find earliest TOI within remaining
            float slice = remaining;
            boolean haveTOI = false;

            if (enableCCD && toiCount < maxTOISteps) {
                float tHit = remaining;
                RigidBody hitA = null;
                RigidBody hitB = null;

                for (RigidBody bullet : bodies) {
                    if (bullet.isStatic() || !bullet.awake) continue;
                    if (!bullet.bullet) continue;

                    // Query nearby candidates using swept AABB (cached per bullet per slice)
                    SweepCacheEntry entry = sweepCache.get(bullet);
                    if (entry == null) {
                        entry = new SweepCacheEntry();
                        sweepCache.put(bullet, entry);
                    }

                    // reuse cache if dt and bullet state unchanged
                    if (!(entry.dt == remaining &&
                          entry.px == bullet.position.x && entry.py == bullet.position.y && entry.pz == bullet.position.z &&
                          entry.vx == bullet.linearVelocity.x && entry.vy == bullet.linearVelocity.y && entry.vz == bullet.linearVelocity.z)) {

                        entry.dt = remaining;
                        entry.px = bullet.position.x; entry.py = bullet.position.y; entry.pz = bullet.position.z;
                        entry.vx = bullet.linearVelocity.x; entry.vy = bullet.linearVelocity.y; entry.vz = bullet.linearVelocity.z;

                        sweptAABBFill(bullet, remaining, entry.aabb);
                        broadphase.queryAABBFill(entry.aabb, entry.candidates, querySeen);
                    }

                    for (RigidBody other : entry.candidates) {
                        if (other == bullet) continue;
                        if (!other.awake && !other.isStatic()) continue;

                        CCD.TOIResult toi = CCD.timeOfImpact(bullet, other, remaining, ccdUseAngularBound);
                        if (toi.hit && toi.t < tHit) {
                            tHit = toi.t;
                            hitA = bullet;
                            hitB = other;
                        }
                    }
                }

                if (hitA != null) {
                    slice = Math.max(tHit, 0f);
                    haveTOI = true;
                    hadAnyTOI = true;
                    toiCount++;

                    if (!hitA.isStatic()) hitA.setAwake(true);
                    if (hitB != null && !hitB.isStatic()) hitB.setAwake(true);
                }
            }

            if (slice < minSubstep) slice = Math.min(minSubstep, remaining);

            // Solve for this slice
            for (int it = 0; it < velocityIterations; it++) {
                for (Joint j : joints) j.solveVelocity(slice);
                for (ContactManifold m : activeManifolds) solveContactVelocity(m, slice);
            }

            // Integrate slice
            integrate(slice);

            // Incremental broadphase updates (moved bodies only)
            flushMovedToBroadphase();

            // Per-substep position correction around TOI
            if (enableCCD && hadAnyTOI && substepPositionIterations > 0) {
                buildContacts(0f);
                for (int it = 0; it < substepPositionIterations; it++) {
                    for (Joint j : joints) j.solvePosition();
                    for (ContactManifold m : activeManifolds) solveContactPosition(m);
                }
                flushMovedToBroadphase();
            }

            remaining -= slice;

            if (!haveTOI) break;

            if (toiCount >= maxTOISteps) {
                if (remaining > 0f) {
                    buildContacts(remaining);
                    for (ContactManifold m : activeManifolds) warmStartContacts(m);
                    for (Joint j : joints) j.warmStart();

                    for (int it = 0; it < velocityIterations; it++) {
                        for (Joint j : joints) j.solveVelocity(remaining);
                        for (ContactManifold m : activeManifolds) solveContactVelocity(m, remaining);
                    }

                    integrate(remaining);
                    flushMovedToBroadphase();
                    remaining = 0f;
                }
                break;
            }
        }

        // position correction (end of frame)
        for (int it = 0; it < positionIterations; it++) {
            for (Joint j : joints) j.solvePosition();
            for (ContactManifold m : activeManifolds) solveContactPosition(m);
        }

        updateSleeping(dt);
    }

    private void markMoved(RigidBody b) {
        if (b == null || b.isStatic() || !b.awake) return;
        if (movedSet.put(b, Boolean.TRUE) == null) movedList.add(b);
    }

    private void flushMovedToBroadphase() {
        if (!incrementalBroadphase) {
            movedSet.clear();
            movedList.clear();
            return;
        }
        // ✅ Ensure cached AABBs are current before updating broadphase
        for (RigidBody b : movedList) {
            b.updateAABB();
            broadphase.update(b);
        }
        movedSet.clear();
        movedList.clear();
    }

    private void integrate(float dt) {
        for (RigidBody b : bodies) {
            if (b.isStatic() || !b.awake) continue;

            // position += v * dt (no allocation)
            b.position.addScaledLocal(b.linearVelocity, dt);

            // ✅ allocation-free quaternion integration:
            // dq = [1, omega*0.5*dt] normalized, then q = dq * q
            float halfDt = 0.5f * dt;
            qDQ.w = 1f;
            qDQ.x = b.angularVelocity.x * halfDt;
            qDQ.y = b.angularVelocity.y * halfDt;
            qDQ.z = b.angularVelocity.z * halfDt;
            qDQ.normalizeLocal(QUAT_NORM_EPS);

            // qNext = dq * q
            qDQ.mulInto(qNext, b.orientation);
            qNext.normalizeLocal(QUAT_NORM_EPS);

            // write back without allocating a new Quat
            b.orientation.set(qNext);

            // Update cached AABB now (so sweptAABB queries and overlap tests see the latest)
            b.updateAABB();

            markMoved(b);
        }
    }

    /** Build contacts using current broadphase grid (no full rebuild). */
    private void buildContacts(float dtHorizon) {
        activeManifolds.clear();

        for (SpatialHashBroadphase3D.Pair p : broadphase.computePairs()) {
            if (!p.a.awake && !p.b.awake) continue;

            // Uses cached AABBs (allocation-free) if you applied the RigidBody cachedAabb patch.
            if (!p.a.aabb().overlaps(p.b.aabb())) continue;

            Collision.detectManifoldInto(p.a, p.b, tmpManifold);
            if (tmpManifold.hit) {
                addManifold(p.a, p.b, tmpManifold);
                continue;
            }

            // Speculative contacts
            if (enableSpeculativeContacts) {
                boolean fastRelevant = p.a.bullet || p.b.bullet;

                scratchRelV.set(
                        p.a.linearVelocity.x - p.b.linearVelocity.x,
                        p.a.linearVelocity.y - p.b.linearVelocity.y,
                        p.a.linearVelocity.z - p.b.linearVelocity.z
                );

                float relSpeed = (float) Math.sqrt(scratchRelV.len2());

                if (fastRelevant && relSpeed > 5f) {
                    GJKDistance.Result dres = GJKDistance.distance(p.a, p.b);
                    if (!dres.intersect) {
                        scratchN.set(dres.normal);
                        if (scratchN.len2() < 1e-10f) scratchN.set(1, 0, 0);
                        else normalizeSelf(scratchN);

                        float closing = -(scratchRelV.dot(scratchN));
                        if (!speculativeRequireClosingSpeed || closing > speculativeClosingSpeedThreshold) {
                            float horizon = Math.max(dtHorizon, 1e-4f);
                            float specDist = relSpeed * horizon * speculativeFactor + speculativeSlop;

                            if (dres.distance <= specDist) {
                                tmpSpecManifold.clear();
                                tmpSpecManifold.hit = true;
                                tmpSpecManifold.normal.set(scratchN);
                                tmpSpecManifold.depth = specDist - dres.distance;

                                tmpSpecManifold.setPointCount(1);
                                tmpSpecManifold.pointsA.get(0).set(dres.pointA);
                                tmpSpecManifold.pointsB.get(0).set(dres.pointB);

                                addManifold(p.a, p.b, tmpSpecManifold);
                            }
                        }
                    }
                }
            }
        }
    }

    private void addManifold(RigidBody a, RigidBody b, Collision.ManifoldResult r) {
        if (!a.isStatic() && !a.awake) a.setAwake(true);
        if (!b.isStatic() && !b.awake) b.setAwake(true);

        long key = SpatialHashBroadphase3D.pairHash(a, b);
        ContactManifold m = manifoldCache.get(key);
        if (m == null) {
            m = new ContactManifold(a, b);
            manifoldCache.put(key, m);
        }

        m.updateFrom(r);
        if (!m.points.isEmpty()) activeManifolds.add(m);
    }

    private void sweptAABBFill(RigidBody b, float dt, AABB out) {
        // Uses cached AABB reference (allocation-free).
        AABB a0 = b.aabb();

        float dx = b.linearVelocity.x * dt;
        float dy = b.linearVelocity.y * dt;
        float dz = b.linearVelocity.z * dt;

        float minx1 = a0.min.x + dx;
        float miny1 = a0.min.y + dy;
        float minz1 = a0.min.z + dz;
        float maxx1 = a0.max.x + dx;
        float maxy1 = a0.max.y + dy;
        float maxz1 = a0.max.z + dz;

        out.min.set(
                Math.min(a0.min.x, minx1),
                Math.min(a0.min.y, miny1),
                Math.min(a0.min.z, minz1)
        );
        out.max.set(
                Math.max(a0.max.x, maxx1),
                Math.max(a0.max.y, maxy1),
                Math.max(a0.max.z, maxz1)
        );

        float s = 0.05f;
        out.min.x -= s; out.min.y -= s; out.min.z -= s;
        out.max.x += s; out.max.y += s; out.max.z += s;
    }

    // --------------------------
    // Allocation-free warm starting
    // --------------------------
    private void warmStartContacts(ContactManifold m) {
        RigidBody A = m.A, B = m.B;
        if ((!A.awake && !A.isStatic()) && (!B.awake && !B.isStatic())) return;

        Vec3 n = m.normal;

        for (ContactPoint cp : m.points) {
            A.worldPointFromLocalInto(sWA, cp.localA);
            B.worldPointFromLocalInto(sWB, cp.localB);

            sP.set(
                    (sWA.x + sWB.x) * 0.5f,
                    (sWA.y + sWB.y) * 0.5f,
                    (sWA.z + sWB.z) * 0.5f
            );
            cp.worldPoint.set(sP);

            sPn.set(n.x * cp.normalImpulse, n.y * cp.normalImpulse, n.z * cp.normalImpulse);
            sPt.set(cp.tangentDir.x * cp.tangentImpulse,
                    cp.tangentDir.y * cp.tangentImpulse,
                    cp.tangentDir.z * cp.tangentImpulse);

            sImp.set(sPn.x + sPt.x, sPn.y + sPt.y, sPn.z + sPt.z);

            sNegImp.set(-sImp.x, -sImp.y, -sImp.z);
            A.applyImpulse(sNegImp, sP);
            B.applyImpulse(sImp, sP);
        }
    }

    // --------------------------
    // Allocation-free contact velocity solver
    // --------------------------
    private void solveContactVelocity(ContactManifold m, float dt) {
        RigidBody A = m.A, B = m.B;
        Vec3 n = m.normal;

        float pen = Math.max(m.penetration - biasSlop, 0f);
        float bias = (pen > 0f) ? (-baumgarte / dt) * pen : 0f;

        for (ContactPoint cp : m.points) {
            Vec3 p = cp.worldPoint;

            sRA.set(p.x - A.position.x, p.y - A.position.y, p.z - A.position.z);
            sRB.set(p.x - B.position.x, p.y - B.position.y, p.z - B.position.z);

            sOmegaXRA.setCross(A.angularVelocity, sRA);
            sVA.set(A.linearVelocity.x + sOmegaXRA.x, A.linearVelocity.y + sOmegaXRA.y, A.linearVelocity.z + sOmegaXRA.z);

            sOmegaXRB.setCross(B.angularVelocity, sRB);
            sVB.set(B.linearVelocity.x + sOmegaXRB.x, B.linearVelocity.y + sOmegaXRB.y, B.linearVelocity.z + sOmegaXRB.z);

            sRV.set(sVB.x - sVA.x, sVB.y - sVA.y, sVB.z - sVA.z);

            float vn = sRV.dot(n);
            float e = (vn < -restitutionThreshold) ? m.restitution : 0f;

            sRaXn.setCross(sRA, n);
            sRbXn.setCross(sRB, n);

            invInertiaWorldMulInto(sInvI_RaXn, A, sRaXn);
            sTmpCross.setCross(sInvI_RaXn, sRA);
            float termA = n.dot(sTmpCross);

            invInertiaWorldMulInto(sInvI_RbXn, B, sRbXn);
            sTmpCross.setCross(sInvI_RbXn, sRB);
            float termB = n.dot(sTmpCross);

            float kN = A.invMass + B.invMass + termA + termB;
            float invKN = 1f / Math.max(1e-8f, kN);

            float jn = -(((1f + e) * vn) + bias) * invKN;

            float oldN = cp.normalImpulse;
            cp.normalImpulse = Math.max(oldN + jn, 0f);
            float dN = cp.normalImpulse - oldN;

            sImp.set(n.x * dN, n.y * dN, n.z * dN);

            sNegImp.set(-sImp.x, -sImp.y, -sImp.z);
            A.applyImpulse(sNegImp, p);
            B.applyImpulse(sImp, p);

            // friction
            sOmegaXRA.setCross(A.angularVelocity, sRA);
            sVA.set(A.linearVelocity.x + sOmegaXRA.x, A.linearVelocity.y + sOmegaXRA.y, A.linearVelocity.z + sOmegaXRA.z);

            sOmegaXRB.setCross(B.angularVelocity, sRB);
            sVB.set(B.linearVelocity.x + sOmegaXRB.x, B.linearVelocity.y + sOmegaXRB.y, B.linearVelocity.z + sOmegaXRB.z);

            sRV.set(sVB.x - sVA.x, sVB.y - sVA.y, sVB.z - sVA.z);

            float rvn = sRV.dot(n);
            sT.set(sRV.x - n.x * rvn, sRV.y - n.y * rvn, sRV.z - n.z * rvn);

            float tl2 = sT.len2();
            if (tl2 < 1e-10f) continue;

            float invTl = 1f / (float) Math.sqrt(tl2);
            sT.mulLocal(invTl);

            sRaXt.setCross(sRA, sT);
            sRbXt.setCross(sRB, sT);

            invInertiaWorldMulInto(sInvI_RaXt, A, sRaXt);
            sTmpCross.setCross(sInvI_RaXt, sRA);
            float termTA = sT.dot(sTmpCross);

            invInertiaWorldMulInto(sInvI_RbXt, B, sRbXt);
            sTmpCross.setCross(sInvI_RbXt, sRB);
            float termTB = sT.dot(sTmpCross);

            float kT = A.invMass + B.invMass + termTA + termTB;
            float invKT = 1f / Math.max(1e-8f, kT);

            float jt = -(sRV.dot(sT)) * invKT;

            float maxF = m.friction * cp.normalImpulse;
            float oldT = cp.tangentImpulse;
            cp.tangentImpulse = MathUtil.clamp(oldT + jt, -maxF, maxF);
            float dT = cp.tangentImpulse - oldT;

            cp.tangentDir.set(sT);

            sImp.set(sT.x * dT, sT.y * dT, sT.z * dT);
            sNegImp.set(-sImp.x, -sImp.y, -sImp.z);

            A.applyImpulse(sNegImp, p);
            B.applyImpulse(sImp, p);
        }
    }

    // --------------------------
    // Allocation-free positional correction
    // --------------------------
    private void solveContactPosition(ContactManifold m) {
        RigidBody A = m.A, B = m.B;
        if (A.isStatic() && B.isStatic()) return;

        float slop = 0.01f;
        float percent = 0.25f;

        float pen = Math.max(m.penetration - slop, 0f);
        if (pen <= 0f) return;

        float invMassSum = A.invMass + B.invMass;
        if (invMassSum <= 0f) return;

        float mag = percent * pen / invMassSum;
        sCorrection.set(m.normal.x * mag, m.normal.y * mag, m.normal.z * mag);

        if (!A.isStatic()) {
            A.position.subScaledLocal(sCorrection, A.invMass);
            A.updateAABB();
            markMoved(A);
        }
        if (!B.isStatic()) {
            B.position.addScaledLocal(sCorrection, B.invMass);
            B.updateAABB();
            markMoved(B);
        }
    }

    // Sleeping / islands (unchanged)
    private void updateSleeping(float dt) {
        UnionFind uf = new UnionFind(bodies.size());
        IdentityHashMap<RigidBody, Integer> idx = new IdentityHashMap<>();
        for (int i = 0; i < bodies.size(); i++) idx.put(bodies.get(i), i);

        for (ContactManifold m : activeManifolds) {
            RigidBody A = m.A, B = m.B;
            if (A.isStatic() || B.isStatic()) continue;
            uf.union(idx.get(A), idx.get(B));
        }

        for (Joint j : joints) {
            RigidBody A = j.bodyA();
            RigidBody B = j.bodyB();
            if (A.isStatic() || B.isStatic()) continue;
            Integer ia = idx.get(A), ib = idx.get(B);
            if (ia != null && ib != null) uf.union(ia, ib);
        }

        class IslandStats {
            boolean canSleep = true;
            float minSleepTime = Float.POSITIVE_INFINITY;
            ArrayList<RigidBody> members = new ArrayList<>();
        }

        HashMap<Integer, IslandStats> islands = new HashMap<>();
        for (int i = 0; i < bodies.size(); i++) {
            RigidBody b = bodies.get(i);
            if (b.isStatic()) continue;
            int root = uf.find(i);
            islands.computeIfAbsent(root, k -> new IslandStats()).members.add(b);
        }

        for (IslandStats st : islands.values()) {
            for (RigidBody b : st.members) {
                if (!b.canSleep) { st.canSleep = false; continue; }
                if (!b.awake) continue;

                float lin2 = b.linearVelocity.len2();
                float ang2 = b.angularVelocity.len2();

                if (lin2 > linearSleepTolerance * linearSleepTolerance ||
                    ang2 > angularSleepTolerance * angularSleepTolerance) {
                    b.sleepTime = 0f;
                    st.canSleep = false;
                } else {
                    b.sleepTime += dt;
                    st.minSleepTime = Math.min(st.minSleepTime, b.sleepTime);
                }
            }

            if (st.canSleep && st.minSleepTime >= timeToSleep) {
                for (RigidBody b : st.members) b.setAwake(false);
            }
        }
    }

    private static final class UnionFind {
        int[] p, r;
        UnionFind(int n) {
            p = new int[n];
            r = new int[n];
            for (int i = 0; i < n; i++) { p[i] = i; r[i] = 0; }
        }
        int find(int x) {
            while (p[x] != x) {
                p[x] = p[p[x]];
                x = p[x];
            }
            return x;
        }
        void union(int a, int b) {
            int ra = find(a), rb = find(b);
            if (ra == rb) return;
            if (r[ra] < r[rb]) p[ra] = rb;
            else if (r[ra] > r[rb]) p[rb] = ra;
            else { p[rb] = ra; r[ra]++; }
        }
    }

    // small helper: normalize Vec3 in-place without allocating
    private static void normalizeSelf(Vec3 v) {
        float len2 = v.x * v.x + v.y * v.y + v.z * v.z;
        if (len2 < 1e-12f) return;
        float inv = 1f / (float) Math.sqrt(len2);
        v.set(v.x * inv, v.y * inv, v.z * inv);
    }

    // invI_world * v (allocation-free): R( invI_body * (R^T v) ), invI_body diagonal
    private void invInertiaWorldMulInto(Vec3 out, RigidBody body, Vec3 vWorld) {
        rotateConjugateInto(sBody, body.orientation, vWorld);
        sScaled.set(
                sBody.x * body.invInertiaBody.m00,
                sBody.y * body.invInertiaBody.m11,
                sBody.z * body.invInertiaBody.m22
        );
        rotateInto(out, body.orientation, sScaled);
    }

    private static void rotateInto(Vec3 out, Quat q, Vec3 v) {
        final float qw = q.w, qx = q.x, qy = q.y, qz = q.z;
        final float vx = v.x, vy = v.y, vz = v.z;

        final float tx = 2f * (qy * vz - qz * vy);
        final float ty = 2f * (qz * vx - qx * vz);
        final float tz = 2f * (qx * vy - qy * vx);

        out.set(
                vx + qw * tx + (qy * tz - qz * ty),
                vy + qw * ty + (qz * tx - qx * tz),
                vz + qw * tz + (qx * ty - qy * tx)
        );
    }

    private static void rotateConjugateInto(Vec3 out, Quat q, Vec3 v) {
        final float qw = q.w, qx = q.x, qy = q.y, qz = q.z;
        final float vx = v.x, vy = v.y, vz = v.z;

        final float tx = 2f * (qy * vz - qz * vy);
        final float ty = 2f * (qz * vx - qx * vz);
        final float tz = 2f * (qx * vy - qy * vx);

        out.set(
                vx - qw * tx + (qy * tz - qz * ty),
                vy - qw * ty + (qz * tx - qx * tz),
                vz - qw * tz + (qx * ty - qy * tx)
        );
    }
}
