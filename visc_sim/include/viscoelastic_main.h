#pragma once

#include <vector>
#include <map>
#include "Eigen/Dense"
#include <cmath>
#include <iostream>

struct Particle
{
    Eigen::Vector2d pos, prev_pos, v;
    int nidx, nx, ny;
    void print()
    {
        std::cout << "Pos " << pos.transpose() << std::endl;
        std::cout << "Pos_prev " << prev_pos.transpose() << std::endl;
        std::cout << "v " << v.transpose() << std::endl;
        std::cout << std::endl;
    }
};

struct Spring
{
    double l; //rest length of spring
    int p_idx1, p_idx2;
};

class FluidSolver
{
    public:
    std::vector<Particle> particles;
    int width, height; //world dims
    
    const double H = 6.4;
    const double RHO0 = 15;
    const double K = 0.5;//15;
    const double KNEAR = 5;//165;
    const double SIGMA = 0;
    const double BETA  = 0.2;//0.3;

    //keep 2 ping ponging grids
    int active_grid; // indicates which grid is currently in use
    std::vector< std::vector<int> > grid[2]; //contains indeces of particles in a given cell (cell indeces are flattened)
    int grid_width, grid_height;

    std::vector< std::vector<int> > nbrs; //array of neighbour indeces

    std::map<int, Spring> springs; //map of a pair of particles to undeformed length
    const int hash_const = 10000; //FIXME: lets assume for now i will not have more than 10000 particles

    //produces unique key (up to 10k particles) given two particle indeces
    int springKey(int particle_idx1, int particle_idx2)
    {
        if (particle_idx1 < particle_idx2)
            return particle_idx1*hash_const+particle_idx2;
        else
        {
            return particle_idx2*hash_const+particle_idx1;
        }
    }

    FluidSolver(int world_width, int world_height)
    {
        width = world_width;
        height = world_height;

        //resize the grid
        grid_width = ceil(world_width/H);
        grid_height = ceil(world_height/H);

        grid[0].resize(grid_height*grid_width);
        grid[1].resize(grid_height*grid_width);
        active_grid = 0;
        //nbrs.resize(grid_height*grid_width);

    }

    void updateParticleIndex(Particle& p)
    {
        p.nx = floor(p.pos(0)/H);
        p.ny = floor(p.pos(1)/H);
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

        // int debug_idx = 1091;
        // particles[debug_idx].print();
        // for (int z = 0; z < nbrs[debug_idx].size(); z++)
        // {
        //     std::cout << nbrs[debug_idx]
        //     printf("NEIGHBOR %d ", nbrs[debug_idx][z]);
        //     particles[nbrs[debug_idx][z]].print();
        // }
        // printf("\n");
        applyViscosity(dt, iter);
        // particles[debug_idx].print();
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].prev_pos = particles[i].pos;
            particles[i].pos += dt*particles[i].v;
        }
        // particles[debug_idx].print();

        // adjustSprings(dt);
        // applySprings(dt);
        doubleDensityRelaxation(dt);

        clearGrid(other_grid);

        // particles[debug_idx].print();

        for (int i = 0; i < particles.size(); i++)
        {
            //update velocity
            Particle& this_p = particles[i];
            this_p.v = (this_p.pos - this_p.prev_pos)/dt;

            //collision resolution TODO better
            Eigen::Vector2d factor1; factor1 << 0,0;//0, 0.9;
            Eigen::Vector2d factor2; factor2 << 0,0;//0.9, 0;
            if (this_p.pos(0)<0)
            {
                this_p.pos(0) = 0.01;
                this_p.v = this_p.v.array() * factor1.array();
            }
            else if (this_p.pos(0) > width)
            {
                this_p.pos(0) = width-0.01;
                this_p.v = this_p.v.array() * factor1.array();
            }

            if (this_p.pos(1)<0)
            {
                this_p.pos(1) = 0.01;
                this_p.v = this_p.v.array() * factor2.array();
            }
            else if (this_p.pos(1) > height)
            {
                this_p.pos(1) = height-0.01;
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

    void adjustSprings(double dt)
    {
        double GAMMA = 0.1;
        double ALPHA = 0.9;

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
                    auto ret = springs.insert(std::pair<int,Spring>(key, s));
                    double rest_len = (ret.first->second).l;
                    double d = GAMMA*rest_len;
                    if (rij_mag > rest_len + d)
                    {
                        rest_len += dt*ALPHA*(rij_mag - rest_len - d);
                    }
                    else if (rij_mag < rest_len - d)
                    {
                        rest_len -= dt*ALPHA*(rest_len-d-rij_mag);
                    }
                }
            }
        }

        std::vector<int> to_delete;
        for (auto it = springs.begin(); it != springs.end(); it++)
        {
            if ((it->second).l > H)
                to_delete.push_back(it->first); 
        }
        
        for (auto it = 0; it < to_delete.size(); it++)
        {
           springs.erase(to_delete[it]);
        }
    }
    
    void applySprings(double dt)
    {
        double K_SPRING = 0.3;
        for (auto it = springs.begin(); it != springs.end(); it++)
        {
            Spring &s = it->second;
            Particle &p1 = particles[s.p_idx1];
            Particle &p2 = particles[s.p_idx2];
            Eigen::Vector2d rij = p1.pos - p2.pos;
            Eigen::Vector2d D = dt*dt*K_SPRING*(1-s.l/H)*(s.l-rij.norm())*rij.normalized();
            
            p1.pos -= 0.5*D;
            p2.pos += 0.5*D;
        }
        // for (int i = 0; i< particles.size(); i++)
        // {
        //     Particle& this_p = particles[i];
        //     for (int j = 0; j < nbrs[i].size(); j++)
        //     {
        //         if (i < nbrs[i][j]) continue;
        //         Particle& nbr_p = particles[nbrs[i][j]];
        //         Eigen::Vector2d rij = nbr_p.pos - this_p.pos;
        //         Eigen::Vector2d rij_norm = rij.normalized();
        //     }
        // }
    }

// break viscoelastic_main.h:217 if iter == 57
    void applyViscosity(double dt, int iter)
    {
        printf("iter %d\n", iter);
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
                    // {
                    //     std::cout << "before " << std::endl;
                    //     assert(0);
                    // }
                    // if (k==711)
                    // {
                    //     std::cout << "before " << std::endl;
                    //     std::cout << "idx " << k << std::endl;
                    //     std::cout << "in visc nbr" << nbr_p.v.transpose() << std::endl;
                    //     std::cout << "in visc this " << this_p.v.transpose() << std::endl;
                    //     //assert(0);
                    // }
                        //     if (k==1666)
                        // {
                        //     std::cout << "q in visc this " << this_p.v.transpose() << std::endl;
                        // }
                    
                    double u = rij_norm.dot(this_p.v - nbr_p.v);
                    if (u > 0)
                    {
                        // if (nbrs[k][j]==1091  )
                        // {
                        //     std::cout << "before in visc nbr " << nbr_p.v.transpose() << std::endl;
                        //     std::cout << "before in visc this " << this_p.v.transpose() << " " << k << std::endl;
                        // }
                        // if (k == 1091)
                        // {
                        //     std::cout << "before in visc this " << this_p.v.transpose() << std::endl;
                        // }

                        Eigen::Vector2d I = dt*(1-q)*(SIGMA*u+BETA*u*u)*rij_norm;
                        this_p.v -= 0.5*I;
                        nbr_p.v += 0.5*I;

                        // if (nbrs[k][j]==1091  )
                        // {
                        //     std::cout << "after in visc nbr " << nbr_p.v.transpose() << std::endl;
                        // }
                        // if (k == 1091)
                        // {
                        //     std::cout << "after in visc this " << this_p.v.transpose() << std::endl;
                        // }
                    }
                    // if (Eigen::isnan(this_p.v.array()).any())
                    // {
                    //     std::cout << "after " << std::endl;
                    //     std::cout << "idx " << k << std::endl;
                    //     std::cout << "in visc nbr" << nbr_p.v.transpose() << std::endl;
                    //     std::cout << "in visc this " << this_p.v.transpose() << std::endl;
                    //     assert(0);
                    // }

                    // if (nbrs[k][j] == 435)
                    // {
                    //     std::cout << "in visc nbr" << nbr_p.v.transpose() << std::endl;
                    //     std::cout << "in visc this " << this_p.v.transpose() << std::endl;
                    //     //std::cout << "in visc nbr " << nbr_p.pos.transpose() << std::endl;
                    //     std::cout << "in visc u dist " << u << " " << dist << std::endl;
                    // }
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
                if (i <= nbr_idx) continue; //TODO should be <=

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