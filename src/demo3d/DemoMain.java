package demo3d;

import physics3d.*;

public class DemoMain {
    public static void main(String[] args){
        World3D world = new World3D(2.0f);
        world.gravity.set(0, -9.8f, 0);
        world.velocityIterations = 24;
        world.positionIterations = 2;
        world.restitutionThreshold = 1.0f;
        world.baumgarte = 0.22f;
        world.enableCCD = true;

        // Ground
        RigidBody ground = new RigidBody(new BoxShape(40, 1, 40), RigidBody.Type.STATIC, 0, new Vec3(0, -2, 0));
        ground.friction = 1.0f;
        world.add(ground);

        // Stack (lets settle; damping helps)
        for (int i=0;i<10;i++){
            RigidBody box = new RigidBody(new BoxShape(0.5f,0.5f,0.5f), RigidBody.Type.DYNAMIC, 2.0f, new Vec3(0, 1f + i*1.01f, 0));
            box.restitution = 0.05f;
            box.friction = 0.9f;
            box.linearDamping = 0.05f;
            box.angularDamping = 0.05f;
            world.add(box);
        }

        // Hinge door with limits + motor + damping (stabilization)
        RigidBody post = new RigidBody(new BoxShape(0.2f, 1.2f, 0.2f), RigidBody.Type.STATIC, 0, new Vec3(-6, 0, 0));
        world.add(post);
        RigidBody door = new RigidBody(new BoxShape(1.4f, 0.1f, 0.8f), RigidBody.Type.DYNAMIC, 3.0f, new Vec3(-4.6f, 0.2f, 0));
        door.friction = 0.7f;
        door.angularDamping = 1.5f; // key for motor stability
        world.add(door);

        HingeJoint hinge = new HingeJoint(post, door, new Vec3(-6, 0.2f, 0), new Vec3(0,1,0));
        hinge.enableLimit = true;
        hinge.lowerAngle = (float)(-Math.PI/3);
        hinge.upperAngle = (float)( Math.PI/3);
        hinge.enableMotor = true;
        hinge.motorSpeed = 1.2f;
        hinge.maxMotorTorque = 12f;
        hinge.softness = 0.05f;
        world.addJoint(hinge);

        // Slider with limits + motor + damping
        RigidBody rail = new RigidBody(new BoxShape(6f, 0.2f, 0.2f), RigidBody.Type.STATIC, 0, new Vec3(0, 0.5f, -6));
        world.add(rail);
        RigidBody slider = new RigidBody(new BoxShape(0.5f,0.5f,0.5f), RigidBody.Type.DYNAMIC, 2.0f, new Vec3(0, 1.5f, -6));
        slider.linearDamping = 0.1f;
        slider.angularDamping = 0.8f;
        world.add(slider);

        SliderJoint sj = new SliderJoint(rail, slider, new Vec3(0, 1.5f, -6), new Vec3(1,0,0));
        sj.enableLimit = true;
        sj.lowerTranslation = -2.5f;
        sj.upperTranslation =  2.5f;
        sj.enableMotor = true;
        sj.motorSpeed = 2.0f;
        sj.maxMotorForce = 40f;
        sj.softness = 0.05f;
        world.addJoint(sj);

        // Fast bullet convex: a small box shot at a thin wall (full convex CCD)
        RigidBody bulletBox = new RigidBody(new BoxShape(0.25f, 0.25f, 0.25f), RigidBody.Type.DYNAMIC, 1.0f, new Vec3(-20, 1.2f, 2));
        bulletBox.linearVelocity.set(60f, 0f, 0f);
        bulletBox.bullet = true;
        bulletBox.angularDamping = 0.5f;
        world.add(bulletBox);

        // Thin wall
        RigidBody wall = new RigidBody(new BoxShape(0.2f, 2.0f, 2.0f), RigidBody.Type.STATIC, 0, new Vec3(8, 1.0f, 2));
        world.add(wall);

        float dt = 1f/60f;
        for (int i=0;i<900;i++){
            world.step(dt);
            if (i % 90 == 0){
                System.out.printf("t=%.2f bulletBox=%s sliderX=%.2f doorAngVel=%s\n", i*dt, bulletBox.position, slider.position.x, door.angularVelocity);
            }
        }
    }
}
