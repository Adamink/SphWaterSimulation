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
#include "rigidbody.h"
#include "constants.h"
#include "kernels.h"

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

        double dx, dy, dz;
        int dh_ratio;  // dh / dx, meaning the ratio of search radius by particle size
        double dh;  // Search radius
        double dt;	// Delta time

        // Including 6 char, each denoting the boundary type of Left/Right/Forward/Back/Down/Up,
        // where 'P' stands for Periodical and 'W' stands for Wall.
        std::string boundary_band;

        int nodes_num;  // The number of all the particles
        // Store liquid nodes, wall nodes and rigidbody nodes together.
        std::vector<Particle> nodes;
        // Store all the index of particles in corresponding grids.
        std::vector<int> idx_table[l + 1][m + 1][h + 1];
        Wheel rigidbody;
        std::vector<Vec3> prev_velocity;
        // Store the number of particles in each grid.
        int grid_num[l + 1][m + 1][h + 1];
        int cur_step;

        // Store vertex coordinates of rigidbody in local coordinate, 
        // where center is the wheel center.
        std::vector<Vec3> localcoord_rigidbody;
        // Store faces of rigidbody with vertex index.
        std::vector<std::vector<int> > faces_rigidbody;
        Vec3 translate_rigidbody;

        Kernels kernels;
    
        Sph();

        void initNodes();
        void setKernels();
        void step();
        void readFileAndShow();

        void compute();
        void computePressureAccelerate(int k);
        void computeVisAccerlerate(int k);
        void computeDensity(int k);
        void updateRigidBody();

        void addNode(int i, int j, int r, std::string status);
        void setWallNodes();
        void setRigidBodyNodes();
        void setLiquidNodes();
        void setRigidBodyMesh();

        void dumpFiles();
        std::string getFilePath(std::string command);
        void dumpLiquidAsCfg(std::string file_path);
        void dumpLiquidAsPly(std::string file_path);
        void dumpRigidBody(std::string file_path);
        bool readLiquidFromCfg(std::string file_path);

        void draw();
        void drawBoundingBox();
        void drawRigidBody();
        void drawLiquid();

        void recordVelocity();
        double getError();
        int getStep();
    };

}
#endif