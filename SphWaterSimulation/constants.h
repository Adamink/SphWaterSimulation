#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace sph_const
{
    const double kXBound = 3.0;
    const double kYBound = 6.0;
    const double kZBOund = 4.0;

    const int kXParticleNum = 30;
    const int kYParticleNum = 60;
    const int kZParticleNum = 40;

    const double kGravity = -9.8;
}

namespace math_const
{
    const double PI = 3.14159265359;
}

namespace camera_const
{
    // Projection parameters.
    const double fovy = 50, aspect = 1, zFar = 1.0 + sph_const::kZParticleNum, zNear = -1.0;
    const double eyex = 3 * sph_const::kXParticleNum;
    const double eyey = 2 * sph_const::kYParticleNum;
    const double eyez = 0.5 * sph_const::kZParticleNum;

    const double centerx = sph_const::kXParticleNum / 2;
    const double centery = sph_const::kYParticleNum / 2;
    const double centerz = sph_const::kZParticleNum / 2;
    const double upx = 0., upy = 0., upz = 1.0;
}

#endif
