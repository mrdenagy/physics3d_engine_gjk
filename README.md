
# Physics3D Java Engine (Near‑Zero GC Build)

A compact but capable **3D rigid‑body physics engine** written in **pure Java** (no external libraries).  
This build is tuned for **near‑zero garbage collection (GC)** during simulation by using **in‑place “Into” APIs**, cached AABBs, pooled broadphase pairs, and scratch workspaces.

---

## Highlights

### Core Simulation
- **3D rigid bodies**: position + quaternion orientation
- **Forces/torques**, gravity, damping
- **Impulse-based solver** with friction + restitution
- **Baumgarte stabilization** + configurable restitution threshold
- **Sleeping + islands** to reduce jitter and improve performance

### Collision Detection
- **Broadphase**: 3D spatial hash with:
  - incremental updates
  - cached candidate pairs
  - reusable long hash set for dedupe (no boxing)
- **Narrowphase**:
  - **Box–Box**: SAT + face clipping → up to 4 contact points (stable stacks)
  - **General convex**: **GJK + EPA** for convex shapes

### Contacts
- **Persistent contact manifolds** (warm-startable)
- **Warm starting**: cached normal/tangent impulses per contact point
- **Contact reduction**: chooses stable contact sets (where applicable)

### Constraints (Joints)
- **Distance joint**
- **Hinge (revolute)** with:
  - motor + limits
  - stabilized axis alignment
- **Slider (prismatic)** with:
  - motor + limits
  - perpendicular constraint solve

### CCD / TOI
- **Continuous Collision Detection** (TOI / conservative advancement)
- Optional **angular motion bound** and substepping behavior (depending on your CCD module version)
- Works best when `body.bullet = true` for fast movers

---

## Near‑Zero GC Design (What’s Different in This Build)

To keep runtime allocations extremely low:
- **Support mapping is allocation‑free**  
  `Shape.supportInto(out, dir, pos, orient)` instead of returning new vectors.
- **AABB computation is allocation‑free**  
  `Shape.computeAABBInto(out, pos, orient)` writes into an existing `AABB`.
- **RigidBody caches its AABB**  
  Each body has `cachedAabb` updated in-place via `updateAABB()`.
- **Broadphase stores bounds objects per body**  
  No `new int[]` bounds per update.
- **Broadphase pair generation uses pooling**  
  Pairs are reused; dedupe uses a primitive long set (no `HashSet<Long>` boxing).