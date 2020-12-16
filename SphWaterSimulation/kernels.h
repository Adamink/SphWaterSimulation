#ifndef KERNELS_H_
#define KERNELS_H_

#include "vec3.h"

class Kernels{
public:
    double search_radius, kpoly6, a_kpoly6, vis_kpoly6;
    Kernels(){};
    Kernels(double search_radius, double kpoly6, double a_kpoly6, double vis_kpoly6):
        search_radius(search_radius), kpoly6(kpoly6), a_kpoly6(a_kpoly6), vis_kpoly6(vis_kpoly6){};
    double Muller03Kernel_Basic(Vec3 origin_point, Vec3 detect_point);
    Vec3 Gradient_Muller03Kernel_Basic(Vec3 origin_point, Vec3 detect_point);
    double Muller03Kernel_Pressure(Vec3 origin_point, Vec3 detect_point);
    Vec3 Gradient_Muller03Kernel_Pressure(Vec3 origin_point, Vec3 detect_point);
    Vec3 Gradient_Muller03Kernel_Vis(Vec3 origin_point, Vec3 detect_point);
};

#endif