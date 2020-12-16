#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <string>
#include "vec3.h"

namespace program_const
{
    // Decide to run simulation or to read from pre-computed files.
    const bool READ_FROM_FILES = false;
    // Number of threads for SPH simulation
    const int NUM_THREADS_COMPUTING = 4;
    const bool IF_VISUALIZE = true;

    const std::string DUMP_FOLDER("../output/final/");
    const std::string CFG_SUBFOLDER = DUMP_FOLDER + "liquid_cfg/";
    const std::string PLY_SUBFOLDER = DUMP_FOLDER + "liquid_ply/";
    const std::string RIGIDBODY_SUBFOLDER = DUMP_FOLDER + "rigidbody_ply/";

    const int kTotalStep = 2000;
    const bool kUseRigidBody = true;

    const int kDumpFileInterval = 5;
    const bool kDumpLiquidAsCfg = false;
    const bool kDumpLiquidAsPly = false;
    const bool kDumpRigidBody = false;
}

namespace sph_const
{

    const double kXBound = 3.0;
    const double kYBound = 6.0;
    const double kZBound = 4.0;

    const int kXParticleNum = 30;
    const int kYParticleNum = 60;
    const int kZParticleNum = 40;

    const double kGravity = 9.8;

    const double kViscosity = 0.02;

    const int kdhRatio = 2;
    const double kLambda = 0.4;
    const double kStiffness = 1000.0;
    const double kMassLiquid = 1.0;
    const double kGamma = 1.0;

    const double kRhoRigidbody = 300.0;
    const double kRigidBodyRadiusOutSize = 0.5;
    const double kRigidBodyRadiusInSize = 0.25;
    const double kRigidBodyHeight = 0.8;
    const Vec3 kRigidBodyCenter = Vec3(kXBound / 2, kYBound / 2, 0.);
    const int kRigidBodyLeafNum = 4;
}

namespace math_const
{
    const double PI = 3.14159265359;
}

namespace visual_const
{
    const int kInitWindowSizeX = 600, kInitWindowSizeY = 600;
    const int kInitWindowPositionX = 10, kInitWindowPositionY = 10;

    // Projection parameters
    const double fovy = 50, aspect = 1, zFar = 1.0 + sph_const::kZParticleNum, zNear = -1.0;
    const double eyex = 3 * sph_const::kXParticleNum;
    const double eyey = 2 * sph_const::kYParticleNum;
    const double eyez = 0.5 * sph_const::kZParticleNum;

    const double centerx = sph_const::kXParticleNum / 2;
    const double centery = sph_const::kYParticleNum / 2;
    const double centerz = sph_const::kZParticleNum / 2;
    const double upx = 0., upy = 0., upz = 1.0;

    const double kParticleRenderingSize = 0.4;
}

#endif
