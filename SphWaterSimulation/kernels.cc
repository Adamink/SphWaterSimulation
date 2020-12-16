#include "kernels.h"

double Kernels::Muller03Kernel_Basic(Vec3 origin_point, Vec3 detect_point){
    double r_square_norm = (detect_point - origin_point).sqrNorm();
    if(search_radius * search_radius > r_square_norm){
        return kpoly6 * pow(search_radius * search_radius - r_square_norm, 3);
    }
    else{
        return 0.0;
    }
}

Vec3 Kernels::Gradient_Muller03Kernel_Basic(Vec3 origin_point, Vec3 detect_point){
    double r_square_norm = (detect_point - origin_point).sqrNorm();
    if(search_radius * search_radius > r_square_norm){
        return (detect_point - origin_point) *
            (6.0 * kpoly6 * pow(search_radius * search_radius - r_square_norm, 2));
    }
    else{
        return Vec3(0.0, 0.0, 0.0);
    }
}

double Kernels::Muller03Kernel_Pressure(Vec3 origin_point, Vec3 detect_point){
    double r_norm = (origin_point - detect_point).norm();
    if(search_radius > r_norm){
        return kpoly6 / 3.0 * pow((search_radius - r_norm), 3);
    }
    else{
        return 0.0;
    }
}

Vec3 Kernels::Gradient_Muller03Kernel_Pressure(Vec3 origin_point, Vec3 detect_point){
    double r_norm = (origin_point - detect_point).norm();
    if(r_norm > search_radius){
        return Vec3(0.0, 0.0, 0.0);
    }
    return (detect_point - origin_point).normalized() *
        (a_kpoly6 * pow(search_radius - r_norm, 2));
}

Vec3 Kernels::Gradient_Muller03Kernel_Vis(Vec3 origin_point, Vec3 detect_point){
    double r_norm = (origin_point - detect_point).norm();
    if(r_norm > search_radius){
        return Vec3(0.0, 0.0, 0.0);
    }
    else{
        return (detect_point - origin_point) * (vis_kpoly6 *
            (-1.5 * r_norm / pow(search_radius, 3) + 2.0 / search_radius / search_radius
                - 0.5 * search_radius / pow(r_norm, 3)));
    }

}