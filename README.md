# physics3d_engine_gjk

A lightweight **3D rigid-body physics / collision** module in Java, centered around **GJK (Gilbert–Johnson–Keerthi)** for convex intersection and **EPA (Expanding Polytope Algorithm)** for penetration depth/normal, plus a practical simulation loop with broadphase, warm starting, joints, and optional CCD.

> **Audience:** developers building simulations, games, or teaching demos that need convex collision detection and a simple rigid-body solver.

---

## Highlights

- **Convex collision detection**
  - **GJK** intersection test (convex shapes)
  - **EPA** for penetration depth and contact normal after GJK hit
  - **GJK distance** for closest points (used by CCD)
  - **Box–Box SAT + face clipping** path producing up to **4 contact points**

- **Broadphase**
  - 3D **spatial hash** broadphase with cached pair generation and allocation-conscious queries

- **Simulation loop** (`World3D`)
  - Semi-implicit integration (forces → velocities)
  - Contact caching via **persistent manifolds** + warm starting
  - Friction + restitution
  - Sleeping / deactivation
  - Optional **CCD** (conservative advancement using GJK distance) for “bullet” bodies
  - Optional speculative contacts (predictive stabilization)

- **Constraints / joints**
  - `HingeJoint` (revolute) with optional **motor + limits**
  - `SliderJoint` (prismatic) with optional **motor + limits**
  - `DistanceJoint`

---

## Project Layout

The source bundle included here is organized as:

```text
main/
  src/
    physics3d/         # engine (math, shapes, collision, broadphase, solver)
    demo3d/            # runnable demo
```

Key entrypoint:
- `demo3d.DemoMain` – builds a small world with stacks, joints, and a fast-moving “bullet” test.

---

## Requirements

- **JDK 17+** recommended (Java 11+ likely works if you avoid newer language features).

---

## Build & Run (no build tool)

From the repository root:

```bash
# compile
mkdir -p out
javac -d out $(find main/src -name "*.java")

# run the demo
java -cp out demo3d.DemoMain
```

Expected output: periodic prints of time, bullet position, slider position, and door angular velocity.

---

## Quick Start (API)

```java
import physics3d.*;

World3D world = new World3D(2.0f);
world.gravity.set(0, -9.8f, 0);
world.enableCCD = true;

RigidBody ground = new RigidBody(
    new BoxShape(40, 1, 40),
    RigidBody.Type.STATIC,
    0,
    new Vec3(0, -2, 0)
);
world.add(ground);

RigidBody box = new RigidBody(
    new BoxShape(0.5f, 0.5f, 0.5f),
    RigidBody.Type.DYNAMIC,
    2.0f,
    new Vec3(0, 2, 0)
);
world.add(box);

float dt = 1f / 60f;
for (int i = 0; i < 600; i++) {
    world.step(dt);
}
```

---

## Docker

This repo snapshot may not include a `Dockerfile` yet. Below is a solid starting point for a **developer-friendly** container build.

### Dockerfile (recommended)

Create `Dockerfile` at repo root:

```dockerfile
# ---- build stage ----
FROM eclipse-temurin:21-jdk AS build
WORKDIR /app

# copy sources
COPY main/src ./main/src

# compile
RUN mkdir -p /app/out \
 && javac -d /app/out $(find main/src -name "*.java")

# ---- run stage ----
FROM eclipse-temurin:21-jre
WORKDIR /app
COPY --from=build /app/out ./out

# default: run demo
CMD ["java", "-cp", "out", "demo3d.DemoMain"]
```

### Build and run

```bash
docker build -t physics3d_engine_gjk:latest .
docker run --rm -it physics3d_engine_gjk:latest
```

### Live iteration (mount your working tree)

```bash
docker run --rm -it \
  -v "$PWD":/app \
  -w /app \
  eclipse-temurin:21-jdk \
  bash -lc 'mkdir -p out && javac -d out $(find main/src -name "*.java") && java -cp out demo3d.DemoMain'
```

### docker-compose (optional)

Create `docker-compose.yml`:

```yaml
services:
  demo:
    build: .
    image: physics3d_engine_gjk:latest
    command: ["java", "-cp", "out", "demo3d.DemoMain"]
```

Run:

```bash
docker compose up --build
```

---

## Notes on Design

### Shapes
- `SphereShape`, `BoxShape`, `ConvexHullShape` implement `Shape`.
- Shapes provide allocation-conscious **support mapping** (`supportInto`) and **AABB** (`computeAABBInto`).

### Narrowphase
- General convex: **GJK → EPA**
- Box–box: **SAT + clipping** to get multiple contact points.

### Broadphase
- `SpatialHashBroadphase3D` hashes AABBs into grid cells, generates unique pairs, and supports query by AABB.

### Solver
- `World3D.step(dt)` integrates forces, builds contacts, warm starts, solves constraints (contacts + joints), and updates sleeping.

---

## Extending

Ideas:
- Add more shapes (capsule, cylinder) using the same `supportInto` pattern.
- Add debug rendering hooks (draw AABBs, contact points, normals).
- Add unit tests for GJK/EPA and joint stability.

---

## License

Add a license file (`LICENSE`) to clarify usage.
