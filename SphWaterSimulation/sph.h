#ifndef SPH_H_
#define SPH_H_

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include <GL/glut.h>

#include "particle.h"
#include "vec3.h"
#include "rigid_body.h"
#include "constants.h"
#include "kernels.h"

using namespace std;

namespace SPH
{
	// Define the size of bounding box of simulation with x in [0, x_bound], y in [0, y_bound], z in [0, z_bound].
	static const double x_bound = sph_const::kXBound, y_bound = sph_const::kYBound, z_bound = sph_const::kZBound;

	// Define the number of particles, and in total l x m x h particles are inited on 3d-grids.
	static const int l = sph_const::kXParticleNum, m = sph_const::kYParticleNum, h = sph_const::kZParticleNum;
	static const Vec3 g_gravity = Vec3(0.0, 0.0, -sph_const::kGravity);

	class Sph{
	public:
		bool use_rigid_body;
		bool if_visualize;
		bool read_from_file;
		int dump_file_interval;

		double max_velocity; // The max velocity of the fluid region
		double mass_liquid;  // The total initial mass of fluid particles
		double rho_liquid;  // The density of liquid
		double init_pressure;  // The initial pressure of particles
		double viscosity;  //  The viscosity of particles

		double lambda, stiffness, KK, gamma;

		int dh_ratio;  // dh / dx, meaning the ratio of search radius by particle size

		// Including 6 char, each denoting the boundary type of Left/Right/Forward/Back/Down/Up,
		// where 'P' stands for Periodical and 'W' stands for Wall.
		std::string boundary_band;

		int nodes_num;  // The number of all the particles
		// Store liquid nodes, wall nodes and rigidbody nodes together.
		vector<Particle> nodes;
		// Store all the index of particles in corresponding grids.
		vector<int> idx_table[l + 1][m + 1][h + 1];
		Wheel rigidbody;
		vector<Vec3> prev_velocity;
		// Store the number of particles in each grid.
		int grid_num[l + 1][m + 1][h + 1];
		int cur_step;

		double dx, dy, dz;
		double dh;  // Search radius
		double dt;	// Delta time

		// Store vertex coordinates of rigidbody in local coordinate, 
		// where center is the wheel center.
		std::vector<Vec3> localcoord_rigidbody;
		// Store faces of rigidbody with vertex index.
		std::vector<vector<int> > faces_rigidbody;
		Vec3 translate_rigidbody;

		Kernels kernels;
		Sph():
			use_rigid_body(program_const::kUseRigidBody),
			if_visualize(program_const::IF_VISUALIZE),
			read_from_file(program_const::READ_FROM_FILES),
			dump_file_interval(program_const::kDumpFileInterval),
			dh_ratio(sph_const::kdhRatio),
			lambda(sph_const::kLambda),
			stiffness(sph_const::kStiffness),
			mass_liquid(sph_const::kMassLiquid),
			gamma(sph_const::kGamma),
			viscosity(sph_const::kViscosity){

			cur_step = 0;
			nodes_num = 0;
			init_pressure = 0.;
			max_velocity = 0.;
			dx = x_bound / l, dy = y_bound / m, dz = z_bound / h;
			dh = dh_ratio * dx;
			rho_liquid = mass_liquid / dx / dy / dz;
			dt = lambda * dh / sqrt(stiffness);
			KK = stiffness * rho_liquid / gamma;
			memset(grid_num, 0, sizeof(grid_num));
			prev_velocity = std::vector<Vec3>(nodes_num, Vec3(0.0, 0.0, 0.0));

			rigidbody = Wheel(sph_const::kRhoRigidbody, sph_const::kRigidBodyRadiusOutSize,
				sph_const::kRigidBodyRadiusInSize, sph_const::kRigidBodyHeight,
				sph_const::kRigidBodyCenter, sph_const::kRigidBodyLeafNum);

			initNodes();
			setKernels();
		}

		void initNodes();
		void setKernels();
		void step();
		void readFileAndShow();

		void compute();
		void computePressureAccelerate(int k);
		void computeVisAccerlerate(int k);
		void computeDensity(int k);

		void setWallNodes();
		void setRigidBodyMesh();
		void setRigidBodyNodes();
		void setLiquidNodes();

		void updateRigidBody();
		void addNode(int i, int j, int r, std::string status);

		std::string getFilePath(std::string command);
		void dumpLiquidAsCfg(std::string file_path);
		void dumpLiquidAsPly(std::string file_path);
		void dumpRigidBody(std::string file_path);
		bool readLiquidFromCfg(std::string file_path);

		void draw();
		void drawBoundingBox();
		void drawRigidBody();
		void drawLiquid();

		double getError();
		void recordVelocity();

		int getStep();

	};

}
#endif