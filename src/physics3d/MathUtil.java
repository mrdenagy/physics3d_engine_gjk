package physics3d;

public final class MathUtil {
    private MathUtil(){}

    public static float clamp(float v, float lo, float hi){
        return Math.max(lo, Math.min(hi, v));
    }

    public static float safeInv(float v){
        return Math.abs(v) < 1e-8f ? 0f : 1f/v;
    }
}
