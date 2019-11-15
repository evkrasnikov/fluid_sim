#pragma once

#include <iostream>
#include "Eigen/Dense"


template <typename T>
struct Vec2d
{
    T x,y;
};

typedef Vec2d<double> Position_t;
typedef Vec2d<double> Velocity_t;
typedef Vec2d<double> Force_t;

//template <typename T> using struct Position Position_t;

// template <typename T>
// struct Velocity
// {
//     T x,y;
// };
// typedef Velocity<double> Velocity_t;

struct Color
{
    uint8_t r,g,b;
};
typedef struct Color Color_t;

//create a struct that holds all the necessary info per particle
//maybe make a template to try different precision values

struct Particle
{
    Eigen::Vector2d x; //position
    Eigen::Vector2d v; //velocity
    Eigen::Vector2d f_pressure; //force
    Eigen::Vector2d f_viscosity; 
    Eigen::Vector2d f_surface;
    Eigen::Vector3i c; //color
    Eigen::Vector2d normal;
    double m; //mass
    double d; //density
    double p; //pressure
    double k; //curvature
};