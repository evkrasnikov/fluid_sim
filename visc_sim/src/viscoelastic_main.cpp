#include "viscoelastic_main.h"
#include <vector>
#include "render.h"
#include <iostream>

#include<SDL2/SDL.h>

#define WORLD_WIDTH 200
#define WORLD_HEIGHT 300

void render_particles(FluidSolver& solver, SDL_Renderer* gRenderer)
{
    for(int k = 0; k < solver.particles.size(); k++)
    {
        Particle& p = solver.particles[k];
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
        SDL_RenderDrawPoint( gRenderer, round(p.pos(0)), round(WORLD_HEIGHT-p.pos(1)) );
    }

}

void init_solver(FluidSolver& solver)
{
    for (int i = 30; i< 70; i++) {
        for (int j = 130; j<200; j++) {
            Particle p;
            p.pos << i,j;
            p.prev_pos << i,j;
            p.v << 0,0;
            p.nidx = 0;
            p.nx = 0;
            p.ny = 0;
            solver.addParticle(p, solver.active_grid);
        }
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

    FluidSolver solver(WORLD_WIDTH, WORLD_HEIGHT);
    double dt = 0.1;
    init_solver(solver);

    //Start up SDL and create window
    if( !init(gWindow, gRenderer, WORLD_WIDTH, WORLD_HEIGHT) )
    {
        printf( "Failed to initialize!\n" );
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
            }

            // render_particles(*g_other, gRenderer);
            // SDL_Delay( 1000 ); //delay in millisecondes

            std::cout << "iteration " << iter << std::endl;
            solver.solve(dt);
            
            SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
            SDL_RenderClear( gRenderer );

            render_particles(solver, gRenderer);

            // SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
            // for (int i = 0; i < 20; i++)
            //     SDL_RenderDrawPoint( gRenderer, x+i+50, y+i+50 );
            
            SDL_RenderPresent(gRenderer);

            SDL_Delay( 10 ); //delay in millisecondes

            iter++;
        
        }   
    }

    //Free resources and close SDL
    close(gWindow, gRenderer);

    return 0;
}