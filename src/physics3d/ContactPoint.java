package physics3d;

/** A single contact point in a persistent manifold with cached impulses for warm starting. */
public final class ContactPoint {
    public Vec3 localA = new Vec3();
    public Vec3 localB = new Vec3();

    public float normalImpulse = 0f;
    public float tangentImpulse = 0f;

    public Vec3 tangentDir = new Vec3(1,0,0);

    public Vec3 worldPoint = new Vec3();

    public void resetImpulses(){
        normalImpulse = 0f;
        tangentImpulse = 0f;
        tangentDir.set(1,0,0);
    }
}
