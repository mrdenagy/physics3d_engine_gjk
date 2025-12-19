package physics3d;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Box-box contact manifold generator using SAT + face clipping.
 * Produces up to 4 contact points.
 *
 * Includes contact reduction heuristics to choose a stable set of 4 points.
 */
public final class BoxBoxManifold {
    private BoxBoxManifold() {}

    public static final class ManifoldResult {
        public boolean hit;
        public Vec3 normal = new Vec3(0,1,0); // from A -> B
        public float penetration;
        public final ArrayList<Vec3> pointsA = new ArrayList<>(4);
        public final ArrayList<Vec3> pointsB = new ArrayList<>(4);
    }

    public static ManifoldResult collide(RigidBody A, RigidBody B, BoxShape aBox, BoxShape bBox){
        ManifoldResult out = new ManifoldResult();

        Vec3[] Ax = boxAxesWorld(A.orientation);
        Vec3[] Bx = boxAxesWorld(B.orientation);

        Vec3 aExt = aBox.halfExtents;
        Vec3 bExt = bBox.halfExtents;

        Vec3 t = B.position.sub(A.position);

        float bestPen = Float.POSITIVE_INFINITY;
        Vec3 bestAxis = new Vec3(1,0,0);
        int bestType = 0; // 0 faceA, 1 faceB, 2 edge-edge
        int bestIndexA = -1;
        int bestIndexB = -1;

        for (int i=0;i<3;i++){
            Vec3 axis = Ax[i];
            float pen = overlapOnAxis(t, axis, Ax, aExt, Bx, bExt);
            if (pen < 0){ out.hit=false; return out; }
            if (pen < bestPen){ bestPen = pen; bestAxis = axis; bestType = 0; bestIndexA = i; }
        }

        for (int i=0;i<3;i++){
            Vec3 axis = Bx[i];
            float pen = overlapOnAxis(t, axis, Ax, aExt, Bx, bExt);
            if (pen < 0){ out.hit=false; return out; }
            if (pen < bestPen){ bestPen = pen; bestAxis = axis; bestType = 1; bestIndexB = i; }
        }

        for (int i=0;i<3;i++){
            for (int j=0;j<3;j++){
                Vec3 axis = Ax[i].cross(Bx[j]);
                float len2 = axis.len2();
                if (len2 < 1e-10f) continue;
                axis.mulLocal(1f/(float)Math.sqrt(len2));
                float pen = overlapOnAxis(t, axis, Ax, aExt, Bx, bExt);
                if (pen < 0){ out.hit=false; return out; }
                float tieBias = 0.002f;
                if (pen + tieBias < bestPen){ bestPen = pen; bestAxis = axis; bestType=2; bestIndexA=i; bestIndexB=j; }
            }
        }

        if (t.dot(bestAxis) < 0) bestAxis = bestAxis.mul(-1f);

        out.hit = true;
        out.normal.set(bestAxis);
        out.penetration = bestPen;

        if (bestType == 0 || bestType == 1){
            buildFaceContact(out, A, B, aBox, bBox, Ax, Bx, bestType, bestIndexA, bestIndexB);
        } else {
            buildEdgeEdgeContact(out, A, B, aBox, bBox, Ax, Bx, bestIndexA, bestIndexB, bestAxis);
        }

        if (out.pointsA.isEmpty()){
            Vec3 p = A.position.add(B.position).mul(0.5f);
            Vec3 pa = p.sub(bestAxis.mul(bestPen*0.5f));
            Vec3 pb = p.add(bestAxis.mul(bestPen*0.5f));
            out.pointsA.add(pa);
            out.pointsB.add(pb);
        }

        // Reduce to max 4 contact points, favoring stability
        reduceTo4(out, out.normal);

        return out;
    }

    private static Vec3[] boxAxesWorld(Quat q){
        return new Vec3[]{
                q.rotate(new Vec3(1,0,0)),
                q.rotate(new Vec3(0,1,0)),
                q.rotate(new Vec3(0,0,1))
        };
    }

    private static float overlapOnAxis(Vec3 t, Vec3 axis, Vec3[] Ax, Vec3 aExt, Vec3[] Bx, Vec3 bExt){
        float ra = aExt.x * Math.abs(axis.dot(Ax[0])) +
                   aExt.y * Math.abs(axis.dot(Ax[1])) +
                   aExt.z * Math.abs(axis.dot(Ax[2]));

        float rb = bExt.x * Math.abs(axis.dot(Bx[0])) +
                   bExt.y * Math.abs(axis.dot(Bx[1])) +
                   bExt.z * Math.abs(axis.dot(Bx[2]));

        float dist = Math.abs(t.dot(axis));
        return (ra + rb) - dist;
    }

    private static void buildFaceContact(ManifoldResult out,
                                         RigidBody A, RigidBody B,
                                         BoxShape aBox, BoxShape bBox,
                                         Vec3[] Ax, Vec3[] Bx,
                                         int bestType, int bestIndexA, int bestIndexB){

        boolean refIsA = (bestType == 0);
        RigidBody Ref = refIsA ? A : B;
        RigidBody Inc = refIsA ? B : A;
        BoxShape refBox = refIsA ? aBox : bBox;
        BoxShape incBox = refIsA ? bBox : aBox;
        Vec3[] Rx = refIsA ? Ax : Bx;
        Vec3[] Ix = refIsA ? Bx : Ax;

        int refAxisIndex = refIsA ? bestIndexA : bestIndexB;

        Vec3 refNormal = Rx[refAxisIndex];
        Vec3 t = Inc.position.sub(Ref.position);
        float sign = (t.dot(refNormal) >= 0) ? 1f : -1f;
        refNormal = refNormal.mul(sign);

        float refExtent = (refAxisIndex==0)? refBox.halfExtents.x : (refAxisIndex==1? refBox.halfExtents.y : refBox.halfExtents.z);
        Vec3 refFaceCenter = Ref.position.add(refNormal.mul(refExtent));

        int incFace = incidentFaceIndex(Ix, refNormal);
        Vec3[] incVerts = faceVerticesWorld(Inc, incBox, incFace);

        int uIndex = (refAxisIndex + 1) % 3;
        int vIndex = (refAxisIndex + 2) % 3;
        Vec3 u = Rx[uIndex];
        Vec3 v = Rx[vIndex];

        float uExt = (uIndex==0)? refBox.halfExtents.x : (uIndex==1? refBox.halfExtents.y : refBox.halfExtents.z);
        float vExt = (vIndex==0)? refBox.halfExtents.x : (vIndex==1? refBox.halfExtents.y : refBox.halfExtents.z);

        Plane[] planes = new Plane[4];
        planes[0] = new Plane(u,  refFaceCenter.add(u.mul(uExt)));
        planes[1] = new Plane(u.mul(-1f), refFaceCenter.add(u.mul(-uExt)));
        planes[2] = new Plane(v,  refFaceCenter.add(v.mul(vExt)));
        planes[3] = new Plane(v.mul(-1f), refFaceCenter.add(v.mul(-vExt)));

        ArrayList<Vec3> poly = new ArrayList<>(4);
        for (Vec3 p : incVerts) poly.add(p);

        for (Plane pl : planes){
            poly = clipPolygon(poly, pl);
            if (poly.isEmpty()) break;
        }

        // Gather candidate contacts
        ArrayList<Candidate> cand = new ArrayList<>();
        for (Vec3 pInc : poly){
            float separation = (pInc.sub(refFaceCenter)).dot(refNormal);
            if (separation <= 0.02f){
                Vec3 pRef = pInc.sub(refNormal.mul(separation));
                Vec3 c = pRef.add(pInc).mul(0.5f);

                // compute 2D coordinates in reference face plane for reduction
                float cu = c.dot(u);
                float cv = c.dot(v);
                cand.add(new Candidate(pRef, pInc, cu, cv, -separation)); // depth positive
            }
        }

        // Reduce candidates to max 4 using spreading heuristic
        ArrayList<Candidate> reduced = reduceCandidates(cand);

        for (Candidate c : reduced){
            if (refIsA){
                out.pointsA.add(c.pRef);
                out.pointsB.add(c.pInc);
            } else {
                out.pointsA.add(c.pInc);
                out.pointsB.add(c.pRef);
            }
        }

        if (!refIsA){
            out.normal.mulLocal(-1f);
        }
    }

    private static void buildEdgeEdgeContact(ManifoldResult out,
                                            RigidBody A, RigidBody B,
                                            BoxShape aBox, BoxShape bBox,
                                            Vec3[] Ax, Vec3[] Bx,
                                            int axisA, int axisB,
                                            Vec3 normal){
        Vec3 pa = supportEdgePoint(A, aBox, Ax, axisA, normal.mul(-1f));
        Vec3 pb = supportEdgePoint(B, bBox, Bx, axisB, normal);
        Vec3 p = pa.add(pb).mul(0.5f);
        out.pointsA.add(p);
        out.pointsB.add(p);
    }

    private static Vec3 supportEdgePoint(RigidBody body, BoxShape box, Vec3[] axes, int edgeAxisIndex, Vec3 dir){
        int i1 = (edgeAxisIndex + 1) % 3;
        int i2 = (edgeAxisIndex + 2) % 3;

        float s1 = (dir.dot(axes[i1]) >= 0) ? 1f : -1f;
        float s2 = (dir.dot(axes[i2]) >= 0) ? 1f : -1f;

        Vec3 local = new Vec3(0,0,0);
        if (i1==0) local.x = box.halfExtents.x * s1;
        if (i1==1) local.y = box.halfExtents.y * s1;
        if (i1==2) local.z = box.halfExtents.z * s1;

        if (i2==0) local.x = box.halfExtents.x * s2;
        if (i2==1) local.y = box.halfExtents.y * s2;
        if (i2==2) local.z = box.halfExtents.z * s2;

        return body.position.add(body.orientation.rotate(local));
    }

    private static int incidentFaceIndex(Vec3[] incAxes, Vec3 refNormal){
        float minDot = Float.POSITIVE_INFINITY;
        int face = 0;
        for (int i=0;i<3;i++){
            float dPos = incAxes[i].dot(refNormal);
            float dNeg = incAxes[i].mul(-1f).dot(refNormal);
            if (dPos < minDot){
                minDot = dPos;
                face = (i==0)?0:(i==1?2:4);
            }
            if (dNeg < minDot){
                minDot = dNeg;
                face = (i==0)?1:(i==1?3:5);
            }
        }
        return face;
    }

    private static Vec3[] faceVerticesWorld(RigidBody body, BoxShape box, int faceIndex){
        float hx=box.halfExtents.x, hy=box.halfExtents.y, hz=box.halfExtents.z;

        Vec3[] local = new Vec3[4];
        switch (faceIndex){
            case 0: // +X
                local[0]=new Vec3( hx, -hy, -hz);
                local[1]=new Vec3( hx, -hy,  hz);
                local[2]=new Vec3( hx,  hy,  hz);
                local[3]=new Vec3( hx,  hy, -hz);
                break;
            case 1: // -X
                local[0]=new Vec3(-hx, -hy,  hz);
                local[1]=new Vec3(-hx, -hy, -hz);
                local[2]=new Vec3(-hx,  hy, -hz);
                local[3]=new Vec3(-hx,  hy,  hz);
                break;
            case 2: // +Y
                local[0]=new Vec3(-hx,  hy, -hz);
                local[1]=new Vec3( hx,  hy, -hz);
                local[2]=new Vec3( hx,  hy,  hz);
                local[3]=new Vec3(-hx,  hy,  hz);
                break;
            case 3: // -Y
                local[0]=new Vec3(-hx, -hy,  hz);
                local[1]=new Vec3( hx, -hy,  hz);
                local[2]=new Vec3( hx, -hy, -hz);
                local[3]=new Vec3(-hx, -hy, -hz);
                break;
            case 4: // +Z
                local[0]=new Vec3(-hx, -hy,  hz);
                local[1]=new Vec3(-hx,  hy,  hz);
                local[2]=new Vec3( hx,  hy,  hz);
                local[3]=new Vec3( hx, -hy,  hz);
                break;
            default: // -Z
                local[0]=new Vec3( hx, -hy, -hz);
                local[1]=new Vec3( hx,  hy, -hz);
                local[2]=new Vec3(-hx,  hy, -hz);
                local[3]=new Vec3(-hx, -hy, -hz);
                break;
        }

        Vec3[] world = new Vec3[4];
        for (int i=0;i<4;i++) world[i] = body.position.add(body.orientation.rotate(local[i]));
        return world;
    }

    private static final class Plane {
        Vec3 n;
        float d;
        Plane(Vec3 n, Vec3 point){
            this.n = n;
            this.d = n.dot(point);
        }
        float distance(Vec3 p){ return n.dot(p) - d; }
    }

    private static ArrayList<Vec3> clipPolygon(ArrayList<Vec3> poly, Plane plane){
        ArrayList<Vec3> out = new ArrayList<>(poly.size());
        if (poly.isEmpty()) return out;

        int count = poly.size();
        Vec3 prev = poly.get(count-1);
        float prevDist = plane.distance(prev);

        for (int i=0;i<count;i++){
            Vec3 curr = poly.get(i);
            float currDist = plane.distance(curr);

            boolean currIn = currDist <= 0f;
            boolean prevIn = prevDist <= 0f;

            if (currIn && !prevIn) out.add(intersect(prev, curr, prevDist, currDist));
            if (currIn) out.add(curr);
            if (!currIn && prevIn) out.add(intersect(prev, curr, prevDist, currDist));

            prev = curr;
            prevDist = currDist;
        }
        return out;
    }

    private static Vec3 intersect(Vec3 a, Vec3 b, float da, float db){
        float t = da / (da - db);
        return a.add(b.sub(a).mul(t));
    }

    // --- Contact reduction helpers ---

    private static final class Candidate {
        Vec3 pRef;
        Vec3 pInc;
        float u, v;
        float depth;
        Candidate(Vec3 pRef, Vec3 pInc, float u, float v, float depth){
            this.pRef = pRef;
            this.pInc = pInc;
            this.u = u;
            this.v = v;
            this.depth = depth;
        }
    }

    private static ArrayList<Candidate> reduceCandidates(ArrayList<Candidate> cand){
        if (cand.size() <= 4) return cand;

        // 1) deepest point
        int i0 = 0;
        for (int i=1;i<cand.size();i++) if (cand.get(i).depth > cand.get(i0).depth) i0 = i;

        // 2) farthest from i0 in plane
        int i1 = farthestInPlane(cand, i0, -1, -1);

        // 3) maximize triangle area
        int i2 = maxAreaWith2(cand, i0, i1);

        // 4) maximize spread: farthest from set (min distance to selected)
        int i3 = farthestFromSet(cand, new int[]{i0,i1,i2});

        ArrayList<Candidate> out = new ArrayList<>(4);
        out.add(cand.get(i0));
        out.add(cand.get(i1));
        out.add(cand.get(i2));
        out.add(cand.get(i3));

        // Sort by angle around centroid in (u,v) for consistent ordering
        float cu=0, cv=0;
        for (Candidate c : out){ cu += c.u; cv += c.v; }
        cu /= out.size(); cv /= out.size();
        final float fcu = cu, fcv = cv;
        Collections.sort(out, (a,b) -> {
            double aa = Math.atan2(a.v - fcv, a.u - fcu);
            double bb = Math.atan2(b.v - fcv, b.u - fcu);
            return Double.compare(aa, bb);
        });
        return out;
    }

    private static int farthestInPlane(ArrayList<Candidate> cand, int i0, int a, int b){
        float best = -1;
        int bestI = 0;
        float u0 = cand.get(i0).u, v0 = cand.get(i0).v;
        for (int i=0;i<cand.size();i++){
            if (i==i0 || i==a || i==b) continue;
            float du = cand.get(i).u - u0;
            float dv = cand.get(i).v - v0;
            float d2 = du*du + dv*dv;
            if (d2 > best){ best = d2; bestI = i; }
        }
        return bestI;
    }

    private static int maxAreaWith2(ArrayList<Candidate> cand, int i0, int i1){
        float u0=cand.get(i0).u, v0=cand.get(i0).v;
        float u1=cand.get(i1).u, v1=cand.get(i1).v;
        float best = -1;
        int bestI = i0;
        for (int i=0;i<cand.size();i++){
            if (i==i0 || i==i1) continue;
            float u2=cand.get(i).u, v2=cand.get(i).v;
            float area2 = Math.abs((u1-u0)*(v2-v0) - (v1-v0)*(u2-u0));
            if (area2 > best){ best = area2; bestI = i; }
        }
        return bestI;
    }

    private static int farthestFromSet(ArrayList<Candidate> cand, int[] sel){
        float best = -1;
        int bestI = sel[0];
        for (int i=0;i<cand.size();i++){
            boolean skip=false;
            for (int s : sel) if (i==s) { skip=true; break; }
            if (skip) continue;

            float minD2 = Float.POSITIVE_INFINITY;
            for (int s : sel){
                float du = cand.get(i).u - cand.get(s).u;
                float dv = cand.get(i).v - cand.get(s).v;
                float d2 = du*du + dv*dv;
                if (d2 < minD2) minD2 = d2;
            }

            if (minD2 > best){ best = minD2; bestI = i; }
        }
        return bestI;
    }

    private static void reduceTo4(ManifoldResult out, Vec3 normal){
        // Safety cap if any path added more than 4.
        int n = out.pointsA.size();
        if (n <= 4) return;

        // Create candidates using midpoint projection to build (u,v) plane basis.
        // Choose arbitrary tangent basis from normal.
        Vec3 t1;
        if (Math.abs(normal.x) < 0.9f) t1 = new Vec3(0, -normal.z, normal.y);
        else t1 = new Vec3(-normal.z, 0, normal.x);
        t1 = t1.normalized();
        Vec3 t2 = normal.cross(t1).normalized();

        ArrayList<Candidate> cand = new ArrayList<>();
        for (int i=0;i<n;i++){
            Vec3 a = out.pointsA.get(i);
            Vec3 b = out.pointsB.get(i);
            Vec3 c = a.add(b).mul(0.5f);
            float u = c.dot(t1);
            float v = c.dot(t2);
            float depth = Math.abs((b.sub(a)).dot(normal));
            cand.add(new Candidate(a,b,u,v,depth));
        }

        ArrayList<Candidate> reduced = reduceCandidates(cand);
        out.pointsA.clear();
        out.pointsB.clear();
        for (Candidate c : reduced){
            out.pointsA.add(c.pRef);
            out.pointsB.add(c.pInc);
        }
    }
}
