#pragma once
#include "Eigen/Dense"

//for density
//returns a scalar
double Wpoly6(Eigen::Vector2d r, double h);


Eigen::Vector2d dWpoly6(Eigen::Vector2d r, double h);


double d2Wpoly6(Eigen::Vector2d r, double h);


//for force due to pressure calculations
//returns a gradient
Eigen::Vector2d dWspiky(Eigen::Vector2d r, double h);


//for viscocity force calculations
//returns a constant
double d2Wvisc(Eigen::Vector2d r, double h);
