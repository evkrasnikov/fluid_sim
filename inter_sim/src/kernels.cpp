#include "kernels.h"

//for density
//returns a scalar
double Wpoly6(Eigen::Vector2d r, double h)
{
    double r_mag = r.norm();
    if (r_mag > h)
        return 0;
    
    double v = 315.0/(64*M_PI*std::pow(h,9))*std::pow(h*h - r_mag*r_mag, 3);
    return v;
}

Eigen::Vector2d dWpoly6(Eigen::Vector2d r, double h)
{
    double r_mag = r.norm();
    if (r_mag > h)
        return Eigen::Vector2d::Zero();
    return -945.0/(32.0*M_PI*std::pow(h,9))*r*std::pow(h*h-r_mag*r_mag,2);
}

double d2Wpoly6(Eigen::Vector2d r, double h)
{
    double r_mag = r.norm();
    if (r_mag > h)
        return 0;
    return -945.0/(32.0*M_PI*std::pow(h,9))*(h*h-r_mag*r_mag)*(3*h*h-7*r_mag*r_mag);
}

//for force due to pressure calculations
//returns a gradient
Eigen::Vector2d dWspiky(Eigen::Vector2d r, double h)
{
    double r_mag = r.norm();
    if (r_mag > h)
        return Eigen::Vector2d::Zero();
    return -45.0/(M_PI*std::pow(h,6))*r.normalized()*std::pow(h-r_mag,2);
}


//for viscocity force calculations
//returns a constant
double d2Wvisc(Eigen::Vector2d r, double h)
{
    double r_mag = r.norm();
    if (r_mag > h)
        return 0;
    
    double v = 45.0/(M_PI*std::pow(h,6))*(h - r_mag);
    return v;
}
