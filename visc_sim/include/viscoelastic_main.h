#pragma once

#include <vector>
#include <map>
#include "Eigen/Dense"
#include <cmath>
#include <iostream>
#include "SDL2/SDL_image.h"

struct Particle
{
    //current and previous position, velocity
    Eigen::Vector2d pos, prev_pos, v;
    
    //grid index
    int nidx, nx, ny;
    void print()
    {
        std::cout << "Pos " << pos.transpose() << std::endl;
        std::cout << "Pos_prev " << prev_pos.transpose() << std::endl;
        std::cout << "v " << v.transpose() << std::endl;
        std::cout << std::endl;
    }

    void reset()
    {
        //pos << 0,0;
        //prev_pos << 0,0;
        //v << 0,0;
        nidx = nx = ny = 0;
    }
};

struct Spring
{
    double l; //rest length of spring
    int p_idx1, p_idx2; //particle indeces on both ends of the spring
};

struct Sphere
{
    Eigen::Vector2d pos; //center of sphere
    double r; //sphere radius
};

/*
 Struct for setting the parameters that control simulation and material properties
*/
struct Settings
{
    // Smoothing kernel radius
    double H;

    // Rest density
    double RHO0;
    
    // "far" density constant
    double K;

    // "near" density constant
    double KNEAR;

    // Viscosity parameters
    double SIGMA;
    double BETA;
    
    // Plasticity and elasticity parameters
    double GAMMA;
    double ALPHA;
    double K_SPRING;

    // Parameters for boundary handling
    double K_STICK;
    double D_STICK;
    double MU;

    // Gravity vector
    Eigen::Vector2d G;

    // Flag that initializes springs with distance between particles
    // instead of H
    bool EN_MOLDING;

    // Flag that enables elasticity and plasticity
    bool EN_SPRINGS;
};

// defines the starting configuration of the water blob
struct WaterRect
{
    double width, height; // width and hight of the rectangular water blob
    double x, y; // position of the top left particle in the rectangle
    double step; // horizontal and vertical separation between the particles
    double v_x, v_y; // initial velocity of the water blob
};

class FluidSolver
{
    public:
    std::vector<Particle> particles;
    int width, height;
    
    // simulation parameters
    double H;
    double RHO0;
    double K;
    double KNEAR;
    double SIGMA;
    double BETA; 
    double GAMMA;
    double ALPHA;
    double K_SPRING; 
    double D_STICK;
    double MU;
    Eigen::Vector2d G;
    bool EN_SPRINGS;
    bool EN_MOLDING;
    double K_STICK;

    //keep 2 ping ponging grids
    int active_grid, other_grid; // indicates which grid is currently in use
    std::vector< std::vector<int> > grid[2]; //contains indeces of particles in a given cell (cell indeces are flattened)
    int grid_width, grid_height;

    std::vector< std::vector<int> > nbrs; //array of neighbour indeces

    std::map<int, Spring> springs; //map of a pair of particles to undeformed length
    const int HASH_CONST = 10000; //FIXME: lets assume for now i will not have more than 10000 particles

    std::vector<Sphere> spheres; // contains all collision spheres in the simulation

    //produces unique key (up to HASH_CONST particles) given two particle indeces
    int springKey(int particle_idx1, int particle_idx2)
    {
        if (particle_idx1 < particle_idx2)
            return particle_idx1*HASH_CONST+particle_idx2;
        else
        {
            return particle_idx2*HASH_CONST+particle_idx1;
        }
    }

    FluidSolver(int world_width, int world_height, Settings sim_settings, std::vector<Sphere> sim_spheres, std::vector<WaterRect> water_blobs)
    {
        width = world_width;
        height = world_height;

        // set the simulation parameters
        H = sim_settings.H;//3;//
        RHO0 = sim_settings.RHO0;//10;//
        K = sim_settings.K;//0.04;//
        KNEAR = sim_settings.KNEAR;//0.1; //
        SIGMA = sim_settings.SIGMA;//1;//
        BETA  = sim_settings.BETA; //0.2;
        GAMMA = sim_settings.GAMMA;
        ALPHA = sim_settings.ALPHA;
        K_SPRING = sim_settings.K_SPRING; //0.3;
        D_STICK = sim_settings.D_STICK;
        MU = sim_settings.MU;
        K_STICK = sim_settings.K_STICK;
        G = sim_settings.G;
        EN_MOLDING = sim_settings.EN_MOLDING;
        EN_SPRINGS = sim_settings.EN_SPRINGS;

        //resize the grid
        grid_width = ceil(world_width/H);
        grid_height = ceil(world_height/H);

        grid[0].resize(grid_height*grid_width);
        grid[1].resize(grid_height*grid_width);
        active_grid = 0;
        other_grid = 1;
        //nbrs.resize(grid_height*grid_width);

        // setup the collision spheres
        spheres = sim_spheres;
        
        // set up the initial water blobs
        for (int k = 0; k < water_blobs.size(); k++)
        {
            WaterRect& w = water_blobs[k];
            for (double i = w.y; i < (w.y + w.height); i += w.step)
            {
                for (double j = w.x; j < (w.x + w.width); j += w.step)
                {
                    Particle p;
                    p.reset();
                    p.pos << j,i;
                    p.prev_pos << j,i;
                    p.v << w.v_x, w.v_y;
                    addParticle(p, active_grid);
                }
            }
        }


    }

    void switchGrid()
    {
        active_grid = other_grid;
        other_grid = active_grid == 0 ? 1 : 0;
    }

    void updateParticleIndex(Particle& p)
    {
        p.nx = floor(p.pos(0)/H);
        p.ny = floor(p.pos(1)/H);
        
        //clip indices to always be in range
        p.nx = std::max(0,std::min(p.nx, grid_width-1));
        p.ny = std::max(0,std::min(p.ny, grid_height-1));

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
        //if (p.nx)
        if (!(p.nx >= 0 && p.nx < grid_width))
        {
            printf("particle id x %d, pos x %f, pos y %f\n", p_idx, p.pos(0), p.pos(1));
            printf("particle id x %d, pos x %f, pos y %f\n", p_idx, p.prev_pos(0), p.prev_pos(1));
            printf("nx %d, ny %d, nidx %d\n", p.nx, p.ny, p.nidx);
            assert(0);
        };
        //particles.push_back(p);
        grid[grid_idx][p.nidx].push_back(p_idx);
    }

    void clearGrid(int idx)
    {
        grid[idx].clear();
        grid[idx].resize(grid_height*grid_width);
    }

    void solve(double dt, int iter)
    {
        findNeighbours();
        std::cout << "active grid " << active_grid << std::endl;
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].v += dt*G;
        }

        applyViscosity(dt, iter);

        clearGrid(other_grid);
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].prev_pos = particles[i].pos;
            particles[i].pos += dt*particles[i].v;

            addParticleGrid(particles[i], other_grid, i);
        }
        switchGrid();
        findNeighbours();

        if (EN_SPRINGS)
        {
            adjustSprings(dt, (iter == 0) && EN_MOLDING);
            applySprings(dt);
        }

        doubleDensityRelaxation(dt);

        clearGrid(other_grid);

        for (int i = 0; i < particles.size(); i++)
        {
            Particle& this_p = particles[i];
            //update velocity
            
            this_p.v = (this_p.pos - this_p.prev_pos)/dt;

            //collision resolution 
            collisionWalls(this_p, dt);
            collisionSpheres(this_p, dt);

            //updateParticleIndex(this_p);
            addParticleGrid(this_p, other_grid, i);

        }
        switchGrid();
        //active_grid = active_grid == 0 ? 1 : 0;
    }

    void doubleDensityRelaxation(double dt)
    {
        for (int i = 0; i< particles.size(); i++)
        {
            //calculate near and far densities
            Particle& this_p = particles[i];
            double rho = 0, rho_near = 0;
            for (int j = 0; j < nbrs[i].size(); j++)
            {
                Particle& nbr_p = particles[nbrs[i][j]];
                Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
                double q = rij.norm()/H;
                if (q < 1)
                {
                    double c = 1-q;
                    rho += c*c;
                    rho_near += c*c*c;
                }
            }
            
            double P = K*(rho-RHO0);
            double P_near = KNEAR*rho_near;
            Eigen::Vector2d dx = Eigen::Vector2d::Zero();
            
            for (int j = 0; j < nbrs[i].size(); j++)
            {
                Particle& nbr_p = particles[nbrs[i][j]];
                Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
                double q = rij.norm()/H;
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

    void collisionWalls(Particle& this_p, double dt)
    {
        Eigen::Vector2d factor1; factor1 << 0,MU;
        Eigen::Vector2d factor2; factor2 << MU,0;
        Eigen::Vector2d sticky_impulse; sticky_impulse << 0,0;
        
        //left wall
        double d_i = fabs(this_p.pos(0));
        if (d_i < D_STICK) 
        {
            Eigen::Vector2d norm_vec; norm_vec << 1,0;
            sticky_impulse = dt*K_STICK*d_i*(1-d_i/D_STICK)*norm_vec;
        }
        if (this_p.pos(0)<0)
        {
            this_p.pos(0) = 0.01;
            this_p.v = this_p.v.array() * factor1.array();
        }

        //right wall
        d_i = fabs(width - this_p.pos(0));
        if (d_i < D_STICK)
        {
            Eigen::Vector2d norm_vec; norm_vec << -1,0;
            sticky_impulse = dt*K_STICK*d_i*(1-d_i/D_STICK)*norm_vec;
        }
        if (this_p.pos(0) > width)
        {
            this_p.pos(0) = width-0.01;
            this_p.v = this_p.v.array() * factor1.array();
        }

        //bottom wall 
        d_i = fabs(this_p.pos(1));
        if (d_i < D_STICK)
        {
            Eigen::Vector2d norm_vec; norm_vec << 0,1;
            sticky_impulse = dt*K_STICK*d_i*(1-d_i/D_STICK)*norm_vec;
        }
        if (this_p.pos(1)<0)
        {
            this_p.pos(1) = 0.01;
            this_p.v = this_p.v.array() * factor2.array();
        }
        
        //top wall
        d_i = fabs(height - this_p.pos(1));
        if (d_i < D_STICK)
        {
            Eigen::Vector2d norm_vec; norm_vec << 0,-1;
            sticky_impulse = dt*K_STICK*d_i*(1-d_i/D_STICK)*norm_vec;
        }
        if (this_p.pos(1) > height)
        {
            this_p.pos(1) = height-0.01;
            this_p.v = this_p.v.array() * factor2.array();
        }

        this_p.v -= sticky_impulse;
    }

    void collisionSpheres(Particle& this_p, double dt)
    {
        
        for (int i = 0; i < spheres.size(); i++)
        {
            Eigen::Vector2d normal = this_p.pos - spheres[i].pos;
            Eigen::Vector2d norm_vec = normal.normalized();
            double dist = normal.norm();

            double d_i = fabs(dist - spheres[i].r);
            Eigen::Vector2d sticky_impulse;
            sticky_impulse << 0, 0;
            if (d_i <= D_STICK)
            {
                sticky_impulse = dt*K_STICK*d_i*(1-d_i/D_STICK)*norm_vec; 
            }
             

            if (dist <= spheres[i].r)
            {
                //split velocity vector into normal and tangential compnents
                Eigen::Vector2d v_norm = this_p.v.dot(norm_vec)*norm_vec;
                this_p.v -= v_norm;
                this_p.v *= MU;

                //extract particle from the sphere
                this_p.pos += (spheres[i].r - dist)*norm_vec;
            }

            this_p.v -= sticky_impulse;
        }
    }

    void adjustSprings(double dt, bool is_first_step)
    {
        //first, look thhrough the nbrs of every particle and 
        // add a spring if distance b/w particles is < H
        for (int i = 0; i< particles.size(); i++)
        {
            Particle& this_p = particles[i];
            for (int j = 0; j < nbrs[i].size(); j++)
            {
                if (i < nbrs[i][j]) continue;
                Particle& nbr_p = particles[nbrs[i][j]];
                Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
                double rij_mag = rij.norm();
                double q = rij_mag/H;
                if (q < 1)
                {
                    //check if spring exists and add if it doesnt
                    int key = springKey(i, nbrs[i][j]);
                    Spring s = {H, i, nbrs[i][j]};
                    if (is_first_step) //"molding"
                        s.l = rij_mag;
                    auto ret = springs.insert(std::pair<int,Spring>(key, s));
                }
            }
        }

        int num_deleted = 0;
        std::vector<int> to_delete;
        //loop through all spring and update rest lengths
        for (auto it = springs.begin(); it != springs.end(); it++)
        {
            Particle& p1 = particles[(it->second).p_idx1];
            Particle& p2 = particles[(it->second).p_idx2];
            double& rest_len = (it->second).l;
            double d = GAMMA*rest_len;
            double rij_mag = (p1.pos - p2.pos).norm();
            
            if (rij_mag > rest_len + d)
            {
                rest_len += dt*ALPHA*(rij_mag - rest_len - d);
            }
            else if (rij_mag < rest_len - d)
            {
                rest_len -= dt*ALPHA*(rest_len-d-rij_mag);
            }

            //check if spring needs to be deleted
            if (rest_len > H)
            {
                to_delete.push_back(it->first); 
                num_deleted++;
            }
            
        }

        printf("Number of deleted springs %d\n", num_deleted);
        for (auto it = 0; it < to_delete.size(); it++)
        {
           springs.erase(to_delete[it]);
        }
    }
    
    void applySprings(double dt)
    {
        
        std::cout << "num springs " << springs.size() << std::endl; 
        std::cout << "num particles " << particles.size() << std::endl; 
        for (auto it = springs.begin(); it != springs.end(); it++)
        {
            Spring &s = it->second;
            Particle &p1 = particles[s.p_idx1]; //i
            Particle &p2 = particles[s.p_idx2]; //j
            Eigen::Vector2d rij = p2.pos - p1.pos;
            Eigen::Vector2d D = dt*dt*K_SPRING*(1-s.l/H)*(s.l-rij.norm())*rij.normalized();
            p1.pos -= 0.5*D;
            p2.pos += 0.5*D;
        }
    }

    void applyViscosity(double dt, int iter)
    {
        //printf("iter %d\n", iter);
        std::cout << nbrs.size() << std::endl;
        for (int k = 0; k < nbrs.size(); k++)
        {
            if (nbrs.size() == 0) continue;
            Particle& this_p = particles[k];
            for (int j = 0; j < nbrs[k].size(); j++)
            {
                if (k < nbrs[k][j]) continue;
                Particle& nbr_p = particles[nbrs[k][j]];
                Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
                Eigen::Vector2d rij_norm = rij.normalized();
                double dist = rij.norm();
                double q = dist/H;
                if (q < 1)
                {
                    // if (Eigen::isnan(this_p.v.array()).any())
                   
                    double u = rij_norm.dot(this_p.v - nbr_p.v);
                    if (u > 0)
                    {
                        Eigen::Vector2d I = dt*(1-q)*(SIGMA*u+BETA*u*u)*rij_norm;
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

        // find the grid cell coordinates where neighbours may reside
        int i_from = p.ny - 1 < 0 ? 0 : p.ny-1;
        int i_to = p.ny + 1 >= grid_height ? p.ny : p.ny+1;
        int j_from = p.nx - 1 < 0 ? 0 : p.nx-1;
        int j_to = p.nx + 1 >= grid_width ? p.nx : p.nx+1;
        
        // accumulate indeces of all particles from neighboring cells into p_nbrs
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
        int nbr_thr = H*H;
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
                if (i <= nbr_idx) continue;

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