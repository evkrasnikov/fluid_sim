#include <iostream>
#include <vector>
#include "main.h"
#include "Eigen/Dense"
#include <functional>
#include <cmath> //for PI
#include "kernels.h"

#include<SDL2/SDL.h>

#include "render.h"

#define WORLD_WIDTH 120
#define WORLD_HEIGHT 140

#define K_PRESSURE 200
#define RHO0_PRESSURE 1
#define MU_VISCOSITY 10.0
#define NORMAL_THR 10.0
#define SIGMA_SURFACE 10.0
#define GRAVITY 9.81



typedef std::function<double(Eigen::Vector2d, double)> SmoothingFuncS; //returns scalar
typedef std::function<Eigen::Vector2d(Eigen::Vector2d, double)> SmoothingFuncV; //returns vector


//TODO: think about quadtrees later
class SpatialGrid
{
    public:
    int width, height;
    double cell_size;

    //grid is setup in such a way that it is safe to access up to +/-1 cell outside of the grid
    std::vector<Particle>** grid; 
    int total_particles;

    SpatialGrid(int _world_width, int _world_height, double _smoothing_width)
    {
        cell_size = _smoothing_width;
        width = ceil((_world_width)/cell_size);
        height = ceil((_world_width)/cell_size);

        grid = new std::vector<Particle>*[height+2];
        for (int i=0; i < height+2; i++)
            grid[i] = new std::vector<Particle>[width+2];
        total_particles=0;
    }

    ~SpatialGrid()
    {
        for (int i=0; i < height+2; i++)
            delete[] grid[i];
        delete[] grid;
    }

    void addParticle(Particle p)
    {
        //find the cell to put it in and append
        int cell_x = floor(p.x(0)/cell_size);
        int cell_y = floor(p.x(1)/cell_size);

        grid[cell_y+1][cell_x+1].push_back(p);
    }

    void clear()
    {
        for (int i = 0; i < height+2; i++)
        {
            for (int j = 0; j < width+2; j++)
            {
                grid[i][j].clear();
            }
        }
    }

    std::vector<Particle>& operator() (int row, int col)
    {
        return grid[row+1][col+1];
    }

};

struct ranges
{
    int i_from, i_to;
    int j_from, j_to;
};

void clip_coordinates(int& i_lo, int& i_hi,
    int& j_lo, int& j_hi,
    int max_i, int max_j)
{
    if (i_lo < 0) i_lo = 0;
    if (i_hi > max_i) i_hi = max_i;
    if (j_lo < 0) j_lo = 0;
    if (j_hi > max_j) j_hi = max_j;
}


void interpDensityNbr(SpatialGrid& g, int i, int j, int k,  SmoothingFuncS f, double h)
{
    Particle& p = g(i,j)[k];
    double density = 0;
    // int ii_from = i - 1;
    // int ii_to = i + 1;
    // int jj_from  = j -1;
    // int jj_to = j + 1;
    // clip_coordinates(ii_from, ii_to, jj_from, jj_to, g.height-1, g.width-1);

    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g(i+ii, j+jj); //g.grid[i+ii][j+jj];
            //for (int kk = 0; kk < p_list.size(); kk++)
            for (Particle& pj : p_list)
            {
                //density += p_list[kk].m*f(p.x - p_list[kk].x);
                density += pj.m*f(p.x - pj.x, h);
            }
        }
    }
    // exclude own contribution (better than if statement on every iteration)
    p.d = density ;//- p.m*f(p.x-p.x, h); 

    //also compute pressure as well
    p.p = K_PRESSURE*(p.d - RHO0_PRESSURE);
}

// void interpDensity(SpatialGrid& g, SmoothingFuncS f, double h)
// {
//     //visit every particle in every cell
//     for (int i = 0; i < g.height; i++)
//     {
//         for (int j = 0; j < g.width; j++)
//         {
//             for (int k = 0; k < g.grid[i][j].size(); k++)
//             {
//                 interpDensityNbr(g, i, j, k, f, h);
//             }
//         }
//     }
//     //f(p.x, 3);
//     std::cout << "interp density\n";
// }

// typedef std::function<void(SpatialGrid& g, int i, int j, int k,  SmoothingFunc f, double h)> InterpFunc;

template<typename InterpFunc>
void interp(SpatialGrid& g, InterpFunc ifunc, SmoothingFuncV f, double h)
{
    for (int i = 0; i < g.height; i++)
    {
        for (int j = 0; j < g.width; j++)
        {
            for (int k = 0; k < g(i,j).size(); k++)
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
            for (int k = 0; k < g(i,j).size(); k++)
            {
                ifunc(g, i, j, k, f, h);
            }
        }
    }
}

void interpForcePressure(SpatialGrid& g, int i, int j, int k,  SmoothingFuncV f, double h)
{
    Particle& p = g(i,j)[k];
    Eigen::Vector2d f_pressure = Eigen::Vector2d::Zero();
    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g(i+ii,j+jj);
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
    Particle& p = g(i,j)[k];
    Eigen::Vector2d f_viscosity = Eigen::Vector2d::Zero();
    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g(i+ii,j+jj);
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
    Particle& p = g(i,j)[k];
    Eigen::Vector2d normal = Eigen::Vector2d::Zero();
    for (int ii = -1; ii <= 1; ii++) //vertical cell offset
    {
        for (int jj = -1; jj <= 1; jj++) //horizontal cell offset
        {
            //get the list of particles in that grid
            std::vector<Particle>& p_list = g(i+ii,j+jj);
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
    Particle& p = g(i,j)[k];
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
            std::vector<Particle>& p_list = g(i+ii,j+jj);
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


//integrate using leapfrog method
void integrate(SpatialGrid& gfrom, SpatialGrid& gto, double dt)
{
    Eigen::Vector2d gravity;
    gravity << 0.0 , -9.81;

    gto.clear();
    for (int i = 0; i < gfrom.height; i++)
    {
        for (int j = 0; j < gfrom.width; j++)
        {
            for(int k = 0; k < gfrom(i,j).size(); k++)
            {
                Particle& p = gfrom(i,j)[k];
                // std::cout << "before: " << p.x.transpose() << std::endl;
                // std::cout << p.v.transpose() << std::endl;
                // std::cout << p.a_prev.transpose() << std::endl;
                Eigen::Vector2d f_tot = p.f_pressure + p.f_surface + p.f_viscosity + gravity*p.m;
                // std::cout << "density " << p.d <<std::endl;
                // std::cout << f_tot.transpose() <<std::endl;
                // std::cout << p.f_pressure.transpose() <<std::endl;
                // std::cout << p.f_surface.transpose() <<std::endl;

                p.x = p.x + p.v*dt +0.5*p.a_prev*dt*dt;
                p.v_prev = p.v;
                p.v = p.v_prev + 0.5*(p.a_prev + f_tot)*dt;
                p.a_prev = f_tot; 

                //detect if particle collided with any surface
                if (p.x(0) < 0)
                {
                    p.x(0) = 0;
                    p.v_prev(0) = -p.v_prev(0);
                    p.v(0) = -p.v(0);
                }
                if (p.x(1) < 0)
                {
                    p.x(1) = 0;
                    p.v_prev(1) = -p.v_prev(1);
                    p.v(1) = -p.v(1);
                    //std::cout << "hi\n";
                }
                if (p.x(0) > WORLD_WIDTH)
                {
                    p.x(0) = WORLD_WIDTH;
                    p.v_prev(0) = -p.v_prev(0);
                    p.v(0) = -p.v(0);
                }
                if (p.x(1) > WORLD_HEIGHT)
                {
                    p.x(1) = WORLD_HEIGHT;
                    p.v_prev(1) = -p.v_prev(1);
                    p.v(1) = -p.v(1);
                }
                // std::cout << "after: " << p.x.transpose() << std::endl;

                gto.addParticle(p);
            }
        }
    }
}


void render_particles(SpatialGrid& g, SDL_Renderer* gRenderer)
{
    for (int i = 0; i < g.height; i++)
    {
        for (int j = 0; j < g.width; j++)
        {
            for(int k = 0; k < g(i,j).size(); k++)
            {
                Particle& p = g(i,j)[k];
                SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
                SDL_RenderDrawPoint( gRenderer, round(p.x(0)), round(WORLD_HEIGHT-p.x(1)) );
            }
        }
    }
}

void init_grid(SpatialGrid& g)
{
    double offset_x = 50, offset_y = 30;
    double step = 0.8;
    for (int i = 0; i < 30; i++)
    for (int j = 0; j < 30; j++)
    {
        Particle p;
        p.x << offset_x + j*step, offset_y + i*step;
        p.v << 0, 0;
        p.m = 1;
        p.v_prev << 0, 0;
        p.a_prev << 0, 0;
        g.addParticle(p);
    }
}


int main(int argc, char* args[])
{
    //The window we'll be rendering to
    SDL_Window* gWindow = NULL;

    //The image we will load and show on the screen
    SDL_Surface* gHelloWorld = NULL;

    //The window renderer
    SDL_Renderer* gRenderer = NULL;

    int x = 0, y = 0;

    const int h = 2;
    SpatialGrid g1(WORLD_WIDTH, WORLD_HEIGHT, h);
    init_grid(g1);
    SpatialGrid g2(WORLD_WIDTH, WORLD_HEIGHT, h);
    double dt = 0.03;

    //Start up SDL and create window
    if( !init(gWindow, gRenderer, WORLD_WIDTH, WORLD_HEIGHT) )
    {
        printf( "Failed to initialize!\n" );
    }
    else
    {
        //Load media
        if( !loadMedia(gHelloWorld) )
        {
            printf( "Failed to load media!\n" );
        }
		else
		{			
			//Main loop flag
			bool quit = false;

			//Event handler
			SDL_Event e;
            int iter = 0;

			//While application is running
			while( !quit )
			{
				//Handle events on queue
				while( SDL_PollEvent( &e ) != 0 )
				{
					//User requests quit
					if( e.type == SDL_QUIT )
					{
						quit = true;
					}
                    // else if(e.type = SDL_KEYDOWN)
                    // {
                    //     switch(e.key.keysym.sym)
                    //     {
                    //         case SDLK_UP:
                    //         y++;
                    //         break;

                    //         case SDLK_DOWN:
                    //         y--;
                    //         break;

                    //         case SDLK_LEFT:
                    //         x--;
                    //         break;

                    //         case SDLK_RIGHT:
                    //         x++;
                    //         break;

                    //         default:
                    //         break;
                    //     }
                    // }
				}

                // render_particles(*g_other, gRenderer);
                // SDL_Delay( 1000 ); //delay in millisecondes

                std::cout << "iteration " << iter << std::endl;
                //compute the forces
                SpatialGrid* g = iter % 2 == 0 ? &g1 : &g2;
                SpatialGrid* g_other = iter % 2 == 0 ? &g2 : &g1;
                interp(*g, interpDensityNbr, Wpoly6, h);
                interp(*g, interpForcePressure, dWspiky, h);
                interp(*g, interpForceViscosity, d2Wvisc, h);
                interp(*g, interpNormal, dWpoly6, h);
                interp(*g, interpForceSurface, d2Wpoly6, h);

                //forward euler to advance the simulation
                integrate(*g, *g_other, dt);
                
                SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
                SDL_RenderClear( gRenderer );

                render_particles(*g_other, gRenderer);

                // SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
                // for (int i = 0; i < 20; i++)
                //     SDL_RenderDrawPoint( gRenderer, x+i+50, y+i+50 );
                
                SDL_RenderPresent(gRenderer);

                SDL_Delay( 10 ); //delay in millisecondes

                iter++;
			}
		}
    }

    //Free resources and close SDL
    close(gWindow, gRenderer);

    return 0;
}