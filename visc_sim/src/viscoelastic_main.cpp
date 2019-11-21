#include "viscoelastic_main.h"
#include <vector>
#include "render.h"
#include <iostream>
#include <map>

#include<SDL2/SDL.h>

#define WORLD_WIDTH 100.0
#define WORLD_HEIGHT 150.0

#define RENDER_WIDTH 400.0
#define RENDER_HEIGHT 600.0

void render_particles(FluidSolver& solver, SDL_Renderer* gRenderer, SDL_Texture* texture)
{
    Eigen::Vector2d c; //shift vector
    c << 0, WORLD_HEIGHT;
    Eigen::Matrix2d s; //scale matrix
    s << RENDER_WIDTH/WORLD_WIDTH, 0, 0, -RENDER_HEIGHT/WORLD_HEIGHT;
    // std::cout << s << std::endl;
    
    for(int k = 0; k < solver.particles.size(); k++)
    {
        Particle& p = solver.particles[k];
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
        Eigen::Vector2d particle_pos = p.pos;
        Eigen::Vector2d world_pos = s*(particle_pos-c);
        //std::cout << "pos " << particle_pos.transpose() << " " << world_pos.transpose() << std::endl;

        if (texture != NULL)
        {
            int x = round(world_pos(0)) - 4;
            int y = round(world_pos(1)) - 4;
            SDL_Rect srcrect = {0,0,9,9};
            SDL_Rect dstrect = {x,y,9,9};
            SDL_RenderCopy(gRenderer, texture, &srcrect, &dstrect);
        }
        else
        {
            SDL_RenderDrawPoint( gRenderer, round(world_pos(0)), round(world_pos(1)) );
        }
        
        
    }

}

void init_solver(FluidSolver& solver)
{
    for (double i = 0; i< 40; i+=1) {
        for (double j = 100; j<130; j+=1) {
            std::cout << i << " " << j << std::endl;
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
    // std::map <int, int> test;
    // test.insert(std::pair<int,int>(4, 1));
    // test.insert(std::pair<int,int>(6, 0));
    // test.insert(std::pair<int,int>(9, 1));
    // test.insert(std::pair<int,int>(1, 1));
    // test.insert(std::pair<int,int>(2, 0));
    // for (auto it = test.begin(); it!=test.end(); it++)
    // {
    //     std::cout << it->first << " " << it->second << std::endl;
    //     if (it->second == 0)
    //         test.erase(it);
    // }
    // return 0;

    //The window we'll be rendering to
    SDL_Window* gWindow = NULL;

    //The image we will load and show on the screen
    SDL_Surface* gHelloWorld = NULL;

    //The window renderer
    SDL_Renderer* gRenderer = NULL;

    //Current displayed texture
    SDL_Texture* gTexture = NULL;   

    FluidSolver solver(WORLD_WIDTH, WORLD_HEIGHT);
    double dt = 1.0/30.0;
    init_solver(solver);

    //Start up SDL and create window
    if( !init(gWindow, gRenderer, RENDER_WIDTH, RENDER_HEIGHT) )
    {
        printf( "Failed to initialize!\n" );
    }
    else
    {
        //try to load the texture
        gTexture = loadTexture("../../assets/dot9.png", gRenderer);
        if (gTexture == NULL)
        {
            printf("Texture was not loaded. Exiting!\n");
            return 0;
        }

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
            solver.solve(dt, iter);
            
            SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
            SDL_RenderClear( gRenderer );

            render_particles(solver, gRenderer, gTexture);

            // SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
            // for (int i = 0; i < 20; i++)
            //     SDL_RenderDrawPoint( gRenderer, x+i+50, y+i+50 );
            // SDL_Rect srcrect = {0,0,5,5};
            // SDL_Rect dstrect = {0,0,5,5};
            // SDL_RenderCopy(gRenderer, gTexture, &srcrect, &dstrect);

            SDL_RenderPresent(gRenderer);

            SDL_Delay( 10 ); //delay in millisecondes

            iter++;
        
        }   
    }

    //Free resources and close SDL
    close(gWindow, gRenderer);

    return 0;
}