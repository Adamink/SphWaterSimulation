#ifndef RIGIDBODY_H_
#define RIGIDBODY_H_

#include "vec3.h"
#include "mat3.h"
#include "particle.h"
#include "constants.h"

class RigidBody{
public:
    Vec3 position;
    Vec3 velocity;
    Vec3 angular_velocity;
    Vec3 linear_momentum, linear_force;
    Vec3 angular_momentum, torque;
    Mat3 rotation;
    // The inverse of inertia matrix.
    Mat3 inertia_inv;

    RigidBody():
        position(Vec3()),
        velocity(Vec3()),
        angular_velocity(Vec3()),
        linear_momentum(Vec3()),
        linear_force(Vec3()),
        angular_momentum(Vec3()),
        torque(Vec3()),
        rotation(Mat3()),
        inertia_inv(Mat3()){}

    virtual void update(double dt) = 0;
};

class Wheel: public RigidBody{
public:
    // Store the index of rigid body particles
    std::vector<int> particle_indexes;

    double m_wheel_radius_outsize;  // The radius of wheel (out circle)
    double m_wheel_radius_insize;  // The radius of wheel (in circle)
    double m_wheel_height;  // The height of wheel
    int m_leafnum;  // The number of leaves in the wheel
    double swirl_velocity;
    double m_Mass;

    Wheel(){
        position = Vec3();
        m_wheel_radius_outsize = 0.0;
        m_wheel_radius_insize = 0.0;
        m_wheel_height = 0.0;
        m_leafnum = 0;
        m_Mass = 1.0;
        particle_indexes.clear();
        inertia_inv = Mat3();

    }

    Wheel(double rho_rigidbody, double radius_outsize, double radius_insize, double height, Vec3 center, int leaf_num){
        position = center;
        m_wheel_radius_outsize = radius_outsize;
        m_wheel_radius_insize = radius_insize;
        m_wheel_height = height;
        m_leafnum = leaf_num;
        m_Mass = rho_rigidbody * height * 4.0 * math_const::PI * \
            (m_wheel_radius_insize * m_wheel_radius_insize + m_wheel_radius_outsize * m_wheel_radius_outsize) / 2.0;
        particle_indexes.clear();
        inertia_inv = Mat3();

    }

    // Semi-implicit update
    virtual void update(double dt){
        swirl_velocity = math_const::PI / 18.0 / dt;
    }
};

class Sphere:public RigidBody{
public:
    std::vector<int> particle_indexes; // Store the index rigid body nodes
    Vec3 m_sphere_center; // The center of sphere
    double m_sphere_radius; // The radius of sphere
    double m_Mass;

    Sphere(){
        m_sphere_radius = 0.0;
        m_Mass = 1.0;
        m_sphere_center = Vec3();
        particle_indexes.clear();
        inertia_inv = Mat3() * (0.4 * m_Mass * m_sphere_radius * m_sphere_radius);

        position = m_sphere_center;
    }

    Sphere(double rho_rigidbody, double radius, Vec3 center){
        m_sphere_radius = radius;
        m_Mass = rho_rigidbody * (4.0 * math_const::PI / 3.0 * pow(radius, 3));
        m_sphere_center = center;
        particle_indexes.clear();
        inertia_inv = Mat3() * (0.4 * m_Mass * m_sphere_radius * m_sphere_radius);

        position = m_sphere_center;
    }

    virtual void update(double dt){ // Semi-implicit
        Mat3 I_inv = rotation.multiply(inertia_inv).multiply(rotation);

        // Linear force & torque has already updated
        linear_momentum += linear_force * dt;
        angular_momentum += torque * dt;

        // Coriolis force
        torque -= angular_velocity.crossProduct(I_inv.multiply(angular_velocity));
        angular_velocity += I_inv.multiply(torque) * dt;

        velocity = linear_momentum / m_Mass;
        position += velocity * dt;
        m_sphere_center = position;

        Mat3 Omega(std::vector<Vec3>{
            Vec3(0.0, -angular_velocity.getz(), angular_velocity.gety()),
                Vec3(angular_velocity.getz(), 0.0, -angular_velocity.getx()),
                Vec3(-angular_velocity.gety(), angular_velocity.getx(), 0.0)});
        rotation += Omega.multiply(rotation) * dt;
    }
};

#endif