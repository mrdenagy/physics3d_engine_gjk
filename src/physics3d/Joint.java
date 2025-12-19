package physics3d;

/**
 * Constraint (joint) between two bodies.
 *
 * This engine uses a sequential impulse solver:
 * - warmStart(): apply cached impulses
 * - solveVelocity(dt): velocity-level constraint enforcement
 * - solvePosition(): optional position-level correction
 */
public interface Joint {
    RigidBody bodyA();
    RigidBody bodyB();

    void warmStart();
    void solveVelocity(float dt);
    void solvePosition();
}
