#include <iostream>
#include <vector>
#include "main.h"
#include "Eigen/Dense"
#include <functional>
#include <cmath> //for PI
#include "kernels.h"

#define GRID_WIDTH 20
#define GRID_HEIGHT 20

#define K_PRESSURE 1.0
#define RHO0_PRESSURE 1.0
#define MU_VISCOSITY 1.0
#define NORMAL_THR 1.0
#define SIGMA_SURFACE 1.0


typedef std::function<double(Eigen::Vector2d, double)> SmoothingFuncS; //returns scalar
typedef std::function<Eigen::Vector2d(Eigen::Vector2d, double)> SmoothingFuncV; //returns vector

//interpDensityNbr(SpatialGrid& g, int i, int j, int k,  SmoothingFunc f, double h)

double test(Eigen::Vector2d r, double h)
{
    std::cout << "opa i like\n";
    return 0.0;
}








//TODO: think about quadtrees later
class SpatialGrid
{
    public:
    int width, height;
    int cell_size;
    std::vector<Particle>** grid; //each cell stores a linked list of particles
    int total_particles;

    SpatialGrid(int _world_width, int _world_height, int _smoothing_width)
    {
        cell_size = _smoothing_width;
        width = (_world_width + _smoothing_width-1)/cell_size;
        height = (_world_width + _smoothing_width-1)/cell_size;

        grid = new std::vector<Particle>*[height];
        for (int i=0; i < height; i++)
            grid[i] = new std::vector<Particle>[width];
        total_particles=0;
    }

    ~SpatialGrid()
    {
        for (int i=0; i < height; i++)
            delete[] grid[i];
        delete[] grid;
    }

    void addParticle(Particle p)
    {
        //find the cell to put it in and append
        int cell_x = floor(p.x(0)/cell_size);
        int cell_y = floor(p.x(1)/cell_size);
        grid[cell_y][cell_x].push_back(p);
    }

};

void interpDensityNbr(SpatialGrid& g, int i, int j, int k,  SmoothingFuncS f, double h)
{
    Particle& p = g.grid[i][j][k];
    double density = 0;
    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g.grid[i+ii][j+jj];
            //for (int kk = 0; kk < p_list.size(); kk++)
            for (Particle& pj : p_list)
            {
                //density += p_list[kk].m*f(p.x - p_list[kk].x);
                density += pj.m*f(p.x - pj.x, h);
            }
        }
    }
    // exclude own contribution (better than if statement on every iteration)
    p.d = density - p.m*f(p.x-p.x, h); 

    //also compute pressure as well
    p.p = K_PRESSURE*(p.d - RHO0_PRESSURE);
}

void interpDensity(SpatialGrid& g, SmoothingFuncS f, double h)
{
    //visit every particle in every cell
    for (int i = 0; i < g.height; i++)
    {
        for (int j = 0; j < g.width; j++)
        {
            for (int k = 0; k < g.grid[i][j].size(); k++)
            {
                interpDensityNbr(g, i, j, k, f, h);
            }
        }
    }
    //f(p.x, 3);
    std::cout << "interp density\n";
}

// typedef std::function<void(SpatialGrid& g, int i, int j, int k,  SmoothingFunc f, double h)> InterpFunc;

template<typename InterpFunc>
void interp(SpatialGrid& g, InterpFunc ifunc, SmoothingFuncV f, double h)
{
    for (int i = 0; i < g.height; i++)
    {
        for (int j = 0; j < g.width; j++)
        {
            for (int k = 0; k < g.grid[i][j].size(); k++)
            {
                ifunc(g, i, j, k, f, h);
            }
        }
    }
}

template<typename InterpFunc>
void interp(SpatialGrid& g, InterpFunc ifunc, SmoothingFuncS f, double h)
{
    for (int i = 0; i < g.height; i++)
    {
        for (int j = 0; j < g.width; j++)
        {
            for (int k = 0; k < g.grid[i][j].size(); k++)
            {
                ifunc(g, i, j, k, f, h);
            }
        }
    }
}

void interpForcePressure(SpatialGrid& g, int i, int j, int k,  SmoothingFuncV f, double h)
{
    Particle& p = g.grid[i][j][k];
    Eigen::Vector2d f_pressure = Eigen::Vector2d::Zero();
    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g.grid[i+ii][j+jj];
            //for (int kk = 0; kk < p_list.size(); kk++)
            for (Particle& pj : p_list)
            {
                //density += p_list[kk].m*f(p.x - p_list[kk].x);
                f_pressure -= pj.m*(pj.p+p.p)/pj.d*f(p.x - pj.x, h);
            }
        }
    }
    // exclude own contribution (better than if statement on every iteration)
    p.f_pressure = 0.5*(f_pressure + p.m*2.0*p.p/p.d*f(p.x-p.x, h)); 
}


void interpForceViscosity(SpatialGrid& g, int i, int j, int k,  SmoothingFuncS f, double h)
{
    Particle& p = g.grid[i][j][k];
    Eigen::Vector2d f_viscosity = Eigen::Vector2d::Zero();
    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g.grid[i+ii][j+jj];
            //for (int kk = 0; kk < p_list.size(); kk++)
            for (Particle& pj : p_list)
            {
                //density += p_list[kk].m*f(p.x - p_list[kk].x);
                f_viscosity += pj.m*(pj.v-p.v)/pj.d*f(p.x - pj.x, h);
            }
        }
    }
    p.f_viscosity = MU_VISCOSITY*f_viscosity; 
}


void interpNormal(SpatialGrid& g, int i, int j, int k, SmoothingFuncV f, double h)
{
    Particle& p = g.grid[i][j][k];
    Eigen::Vector2d normal = Eigen::Vector2d::Zero();
    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g.grid[i+ii][j+jj];
            //for (int kk = 0; kk < p_list.size(); kk++)
            for (Particle& pj : p_list)
            {
                //density += p_list[kk].m*f(p.x - p_list[kk].x);
                normal += pj.m/pj.d*f(p.x - pj.x, h);
            }
        }
    }
    p.normal = normal;
}

void interpForceSurface(SpatialGrid& g, int i, int j, int k, SmoothingFuncS f, double h)
{
    Particle& p = g.grid[i][j][k];
    double curvature = 0;
    double n_mag = p.normal.norm();
    if (n_mag < NORMAL_THR)
    {
        p.f_surface = Eigen::Vector2d::Zero();
        p.k = 0;
        return;
    }

    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g.grid[i+ii][j+jj];
            //for (int kk = 0; kk < p_list.size(); kk++)
            for (Particle& pj : p_list)
            {
                //density += p_list[kk].m*f(p.x - p_list[kk].x);
                curvature -= pj.m/pj.d*f(p.x - pj.x, h);
            }
        }
    }
    p.k = curvature/n_mag;
    p.f_surface = SIGMA_SURFACE*p.normal*p.k; 
}

//need to create the grid
//the grid must know about the size of the render area
//it must split the render area according to the cell size,
//which should be related to the width parameter of the smoothing kernel
//every cell should contain all particles in that cell (copy or reference?)
//particles should be kept in a linked list or a vector?

//need to keep a datastructure that stores all particles 
//at the beginning of each simulation step need to move particles of this data structure to the grid
//may potentially want to have to grid structure and ping pong between the two

//what to do???
//potentiall create a generic interpolation function and then use it to interp all other quantities
//write the main loop of the simulation with 2 ping ponging grids

void integrate(SpatialGrid& gfrom, SpatialGrid& gto, double dt)
{
    
}

int main(void)
{
    std::cout << "waow!" << std::endl;
    Particle p;
    //interpDensity(p, test);
    //try to create a simulation loop
    const int h = 4;
    SpatialGrid g1(1280, 720, h);
    SpatialGrid g2(1280, 720, h);
    for (int i = 0; i < 100; i++)
    {
        //compute the forces
        SpatialGrid* g = i % 2 == 0 ? &g1 : &g2;
        interp(*g, interpDensityNbr, Wpoly6, h);
        interp(*g, interpForcePressure, dWspiky, h);
        interp(*g, interpForceViscosity, d2Wvisc, h);
        interp(*g, interpNormal, dWpoly6, h);
        interp(*g, interpForceSurface, d2Wpoly6, h);

        //forward euler to advance the simulation


    }
    return 0;
}