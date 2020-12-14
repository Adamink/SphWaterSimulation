#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <string>
#include "vec3.h"

class Particle{
public:
	Vec3 velocity;
	Vec3 position;
	double rho;  // Density of the particle
	double pressure;  // Pressure
	Vec3 force_vis;
	Vec3 force_press;
	Vec3 acc_ext;
	// The status of particle, "Wall" for Wall, "RigidBody" for RigidBody, "Liquid" for Liquid
	std::string status;
	vector<int> neighbor_index;

	Particle(double rho0, double init_pressure, Vec3 acc_g):
		velocity(Vec3(0.0, 0.0, 0.0)),
		position(Vec3(0.0, 0.0, 0.0)),
		force_vis(Vec3(0.0, 0.0, 0.0)),
		force_press(Vec3(0.0, 0.0, 0.0)),
		acc_ext(acc_g),
		rho(rho0),
		pressure(init_pressure),
		status("L"){}

	void setPosition(double new_x, double new_y, double new_z){
		position = Vec3(new_x, new_y, new_z);
	}

	void setPosition(Vec3 new_place){
		position = new_place;
	}

	void setStatus(std::string new_status){
		status = new_status;
	}

	bool isWall(){
		return status == "Wall";
	}
	bool isRigidBody(){
		return status == "RigidBody";
	}
	bool isFluid(){
		return status == "Liquid";
	}
};


#endif