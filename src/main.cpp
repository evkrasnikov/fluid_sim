#include <iostream>
#include <vector>
#include "main.h"
#include "Eigen/Dense"
#include <functional>
#include <cmath> //for PI

#define GRID_WIDTH 20
#define GRID_HEIGHT 20



typedef std::function<double(Eigen::Vector2d, double)> SmoothingFunc;

double test(Eigen::Vector2d r, double h)
{
    std::cout << "opa i like\n";
    return 0.0;
}

double Wpoly6(Eigen::Vector2d r, double h)
{
    double r_mag = r.norm();
    if (r_mag > h)
        return 0;
    
    double v = 315.0/(64*M_PI*std::pow(h,9))*std::pow(h*h - r_mag*r_mag, 3);
    return v;
}



void interpForcePressure(Particle& p)
{

}

void interpForceForceViscosity(Particle& p)
{

}

void interpForceSurface(Particle& p)
{

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

void interpDensityNbr(SpatialGrid& g, int i, int j, int k, Particle p, SmoothingFunc f, double h)
{
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
}

void interpDensity(SpatialGrid& g, SmoothingFunc f, double h)
{
    //visit every particle in every cell
    for (int i = 0; i < g.height; i++)
    {
        for (int j = 0; j < g.width; j++)
        {
            for (int k = 0; k < g.grid[i][j].size(); k++)
            {
                interpDensityNbr(g, i, j, k, g.grid[i][j][k], f, h);
            }
        }
    }
    //f(p.x, 3);
    std::cout << "interp density\n";
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


int main(void)
{
    std::cout << "waow!" << std::endl;
    Particle p;
    //interpDensity(p, test);
    return 0;
}