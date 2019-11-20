#pragma once

#include <vector>
#include "Eigen/Dense"
#include <cmath>
#include <iostream>

struct Particle
{
    Eigen::Vector2d pos, prev_pos, v;
    int nidx, nx, ny;
};

class FluidSolver
{
    public:
    std::vector<Particle> particles;
    int width, height; //world dims
    
    const double h = 3.4;
    const double rho0 = 15;
    const double k = 0.5;
    const double knear = 5;
    const double sigma = 0;
    const double beta  = 0.3;

    //keep 2 ping ponging grids
    int active_grid; // indicates which grid is currently in use
    std::vector< std::vector<int> > grid[2]; //contains indeces of particles in a given cell (cell indeces are flattened)
    int grid_width, grid_height;

    std::vector< std::vector<int> > nbrs; //array of neighbour indeces

    FluidSolver(int world_width, int world_height)
    {
        width = world_width;
        height = world_height;

        //resize the grid
        grid_width = ceil(world_width/h);
        grid_height = ceil(world_height/h);

        grid[0].resize(grid_height*grid_width);
        grid[1].resize(grid_height*grid_width);
        active_grid = 0;
        //nbrs.resize(grid_height*grid_width);

    }

    void updateParticleIndex(Particle& p)
    {
        p.nx = floor(p.pos(0)/h);
        p.ny = floor(p.pos(1)/h);
        p.nidx = p.ny * grid_width + p.nx;;
    }
    void addParticle(Particle& p, int grid_idx)
    {
        updateParticleIndex(p);
        particles.push_back(p);
        grid[grid_idx][p.nidx].push_back(particles.size()-1);
    }
     void addParticleGrid(Particle& p, int grid_idx, int p_idx)
    {
        updateParticleIndex(p);
        //particles.push_back(p);
        grid[grid_idx][p.nidx].push_back(p_idx); //FIXME DUMBASS
    }

    void clearGrid(int idx)
    {
        grid[idx].clear();
        grid[idx].resize(grid_height*grid_width);
    }

    void solve(double dt)
    {
        int other_grid = active_grid == 0 ? 1 : 0;
        Eigen::Vector2d g;
        g << 0, -9.81;
        
        std::cout << "active grid " << active_grid << std::endl;
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].v += dt*g;
        }

        findNeighbours();
        // for (int e = 0; e < nbrs[50].size(); e++)
        //     std::cout << nbrs[50][e] << " ";
        // std::cout << std::endl;
        // std::cout << particles[0].pos << " " << particles[0].nx << " " << particles[0].ny << std::endl;
        // std::cout << nbrs.size() << std::endl;
        // std::cout << "in cell ";
        // for (int e = 0; e < grid[active_grid][7*grid_width+8].size(); e++)
        // {
        //     std::vector<int>& p = grid[active_grid][7*grid_width+8];
        //     std::cout << p[e] << " ";
        // }

        applyViscosity(dt);

        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].prev_pos = particles[i].pos;
            particles[i].pos += dt*particles[i].v;
        }

        doubleDensityRelaxation(dt);

        clearGrid(other_grid);

        for (int i = 0; i < particles.size(); i++)
        {
            //update velocity
            Particle& this_p = particles[i];
            this_p.v = (this_p.pos - this_p.prev_pos)/dt;

            //collision resolution TODO better
            Eigen::Vector2d factor1; factor1 << 0, 0.9;
            Eigen::Vector2d factor2; factor2 << 0.9, 0;
            if (this_p.pos(0)<0)
            {
                this_p.pos(0) = 0;
                this_p.v = this_p.v.array() * factor1.array();
            }
            else if (this_p.pos(0) > width)
            {
                this_p.pos(0) = width;
                this_p.v = this_p.v.array() * factor1.array();
            }

            if (this_p.pos(1)<0)
            {
                this_p.pos(1) = 0;
                this_p.v = this_p.v.array() * factor2.array();
            }
            else if (this_p.pos(1) > height)
            {
                this_p.pos(1) = height;
                this_p.v = this_p.v.array() * factor2.array();
            }

            //updateParticleIndex(this_p);
            addParticleGrid(this_p, other_grid, i);

        }
        active_grid = active_grid == 0 ? 1 : 0;
    }

    void doubleDensityRelaxation(double dt)
    {
        for (int i = 0; i< particles.size(); i++)
        {
            Particle& this_p = particles[i];
            double rho = 0, rho_near = 0;
            for (int j = 0; j < nbrs[i].size(); j++)
            {
                Particle& nbr_p = particles[nbrs[i][j]];
                Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
                double q = rij.norm()/h;
                if (q < 1)
                {
                    double c = 1-q;
                    rho += c*c;
                    rho_near += c*c*c;
                }
            }
            
            double P = k*(rho-rho0);
            double P_near = knear*rho_near;
            Eigen::Vector2d dx = Eigen::Vector2d::Zero();
            
            for (int j = 0; j < nbrs[i].size(); j++)
            {
                Particle& nbr_p = particles[nbrs[i][j]];
                Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
                double q = rij.norm()/h;
                if (q < 1)
                {
                    double c = 1-q;
                    Eigen::Vector2d D = dt*dt*(P*c+P_near*c*c)*rij.normalized();
                    nbr_p.pos += 0.5*D;
                    dx -= 0.5*D;
                }
            }
            this_p.pos += dx;
        }
        
    }

    void applyViscosity(double dt)
    {
        for (int k = 0; k < nbrs.size(); k++)
        {
            if (nbrs.size() == 0) continue;
            Particle& this_p = particles[k];
            for (int j = 0; j < nbrs[k].size(); j++)
            {
                Particle& nbr_p = particles[nbrs[k][j]];
                Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
                Eigen::Vector2d rij_norm = rij.normalized();
                double dist = rij.norm();
                double q = dist/h;
                if (q < 1)
                {
                    double u = rij_norm.dot(this_p.v - nbr_p.v);
                    if (u > 0)
                    {
                        Eigen::Vector2d I = dt*(1-q)*(sigma*u+beta*u*u)*rij_norm;
                        this_p.v -= 0.5*I;
                        nbr_p.v += 0.5*I;
                    }
                }
            }
        }
    }

    void getPotentialNbrs(std::vector<int>& p_nbrs, Particle& p)
    {
        p_nbrs.clear();
        int i_from = p.ny - 1 < 0 ? 0 : p.ny-1;
        int i_to = p.ny + 1 >= grid_height ? p.ny : p.ny+1;
        int j_from = p.nx - 1 < 0 ? 0 : p.nx-1;
        int j_to = p.nx + 1 >= grid_width ? p.nx : p.nx+1;
        for (int ii = i_from; ii <= i_to; ii++)
        for (int jj = j_from; jj <= j_to; jj++)
        {
            int nidx = ii*grid_width+jj;
            p_nbrs.insert(p_nbrs.end(), grid[active_grid][nidx].begin(), grid[active_grid][nidx].end());
        }
    }

    void findNeighbours()
    {
        int num_particles = particles.size();
        int nbr_thr = h*h;
        nbrs.clear();
        nbrs.resize(num_particles);

        for (int i = 0; i < num_particles; i++)
        {
            Particle p = particles[i];
            std::vector<int> p_nbrs;
            getPotentialNbrs(p_nbrs, p);

            for (int j = 0; j < p_nbrs.size(); j++)
            {
                int nbr_idx = p_nbrs[j];
                if (i < nbr_idx) continue;

                Particle p_other = particles[nbr_idx];

                if ((p.pos - p_other.pos).norm() <= nbr_thr)
                {
                    nbrs[i].push_back(nbr_idx);
                    nbrs[nbr_idx].push_back(i);
                }
            }
        }
    }
};