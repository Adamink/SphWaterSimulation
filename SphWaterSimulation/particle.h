#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "vec3.h"

class Particle{
public:
	Vec3 velocity; //粒子速度 | velocity of particle
	Vec3 position; //粒子位置 | position of particle
	double rho; //粒子密度 | density of particle
	double p; //压强 | Pressure
	Vec3 force_vis;
	Vec3 force_press;
	Vec3 acc_ext;
	char status; //粒子的状态属性（固体/液体/气体）| The status of particle (solid/liquid/gas)
	vector<int> neighbor_index;

	Particle(double rho0, double p0, Vec3 acc_g):
		velocity(Vec3(0.0, 0.0, 0.0)),
		position(Vec3(0.0, 0.0, 0.0)),
		force_vis(Vec3(0.0, 0.0, 0.0)),
		force_press(Vec3(0.0, 0.0, 0.0)),
		acc_ext(acc_g),
		rho(rho0),
		p(p0),
		status('L'){}

	//改变粒子的位置信息
	void setplace(double new_x, double new_y, double new_z){
		position = Vec3(new_x, new_y, new_z);
	}

	void setplace(Vec3 new_place){
		position = new_place;
	}

	void set_status(char ch){
		status = ch;
	}

	bool isWall(){
		return status == 'W';
	}
	bool isRigidBody(){
		return status == 'R';
	}
	bool isFluid(){
		return status == 'L';
	}
};


#endif