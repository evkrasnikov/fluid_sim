#include "viscoelastic_main.h"
#include <vector>
#include "render.h"
#include <iostream>
#include <map>

#include<SDL2/SDL.h>
#include<SDL2/SDL_image.h>

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

void render_spheres(FluidSolver& solver, SDL_Renderer* gRenderer, SDL_Texture* texture)
{
    Eigen::Vector2d c; //shift vector
    c << 0, WORLD_HEIGHT;
    Eigen::Matrix2d s_mat; //scale matrix
    s_mat << RENDER_WIDTH/WORLD_WIDTH, 0, 0, -RENDER_HEIGHT/WORLD_HEIGHT;
    // std::cout << s << std::endl;
    for(int k = 0; k < solver.spheres.size(); k++)
    {
        Sphere s = solver.spheres[k];
        //SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
        Eigen::Vector2d sphere_pos = s.pos;
        Eigen::Vector2d world_pos = s_mat*(sphere_pos-c);
        Eigen::Vector2d world_r; 
        world_r << s.r*RENDER_WIDTH/WORLD_WIDTH, s.r*RENDER_HEIGHT/WORLD_HEIGHT;
        //std::cout << "pos " << particle_pos.transpose() << " " << world_pos.transpose() << std::endl;

        int x = round(world_pos(0)-world_r(0));
        int y = round(world_pos(1)-world_r(1));
        SDL_Rect srcrect = {0,0,51,51};
        SDL_Rect dstrect = {x,y,int(world_r(0)*2),int(world_r(1)*2)};
        SDL_RenderCopy(gRenderer, texture, &srcrect, &dstrect);

        
        
    }
}

void set_settings(Settings& set)
{
    set.H = 3.4;//3;//
    set.RHO0 = 15;//10;//
    set.K = 0.5;//0.04;//
    set.KNEAR = 5;//0.1; //
    set.SIGMA = 0;//1;//
    set.BETA  = 0.2; //0.2;
    
    // plasticity and elasticity 
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 0.3;//40; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;

    // boundary stuff
    set.K_STICK = 35;//35;
}

void set_spheres(std::vector<Sphere>& spheres)
{
    Sphere s;
    s.pos << 50,50;
    s.r = 20;
    spheres.push_back(s);
    s.pos << 20, 20;
    s.r = 10;
    spheres.push_back(s);
}

void init_solver(FluidSolver& solver)
{
    for (double i = 70; i< 100; i+=1) {
        for (double j = 30; j<60; j+=1) {
            std::cout << i << " " << j << std::endl;
            Particle p;
            p.pos << j,i;
            p.prev_pos << j,i;
            p.v << 0,0;
            p.nidx = 0;
            p.nx = 0;
            p.ny = 0;
            solver.addParticle(p, solver.active_grid);
        }
    }
}

void init_solver2(FluidSolver& solver)
{
    for (double i = 1; i< 20; i+=1) {
        for (double j = 30; j<70; j+=1) {
            std::cout << i << " " << j << std::endl;
            Particle p;
            p.pos << j,i;
            p.prev_pos << j,i;
            p.v << 0,0;
            p.nidx = 0;
            p.nx = 0;
            p.ny = 0;
            solver.addParticle(p, solver.active_grid);
        }
    }
}

void init_solver3(FluidSolver& solver)
{
    for (double i = 120; i< 140; i+=1) {
        for (double j = 0; j<100; j+=1) {
            Particle p;
            p.pos << j,i;
            p.prev_pos << j,i;
            p.v << 0,0;
            p.nidx = 0;
            p.nx = 0;
            p.ny = 0;
            solver.addParticle(p, solver.active_grid);
        }
    }
}

void init_solver4(FluidSolver& solver)
{
    for (double i = 110; i< 140; i+=1) {
        for (double j = 0; j<50; j+=1) {
            Particle p;
            p.pos << j,i;
            p.prev_pos << j,i;
            p.v << 15,0;
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
    
    //
    //create a directory to store pngs
    system("rm -r ./pngs");
    system("mkdir -p ./pngs");
    char filename_buf[100];

    //The window we'll be rendering to
    SDL_Window* gWindow = NULL;

    //The image we will load and show on the screen
    SDL_Surface* gHelloWorld = NULL;

    //The window renderer
    SDL_Renderer* gRenderer = NULL;

    //Current displayed texture
    SDL_Texture* textureParticle = NULL;
    SDL_Texture* textureSphere = NULL;   

    SDL_Surface *frm = SDL_CreateRGBSurface(0, RENDER_WIDTH, RENDER_HEIGHT, 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000); //holds the output frame

    //WaterRect water_blob = {100, 20, 0, 120, 1, 0, 0};
    WaterRect water_blob = {30, 30, 30, 120, 1, -10, 0};
    std::vector<Sphere> spheres;
    set_spheres(spheres);
    Settings set;
    set_settings(set);

    FluidSolver solver(WORLD_WIDTH, WORLD_HEIGHT, set, spheres, water_blob);
    double dt = 1.0/30.0;///4.0;
    //init_solver(solver);
    //init_solver3(solver);

    // struct alol {int x; double y;} a;
    // double &gg = a.y;
    // gg = 6;
    // std::cout << gg << " " << a.y << std::endl;
    // return 0;

    //Start up SDL and create window
    if( !init(gWindow, gRenderer, RENDER_WIDTH, RENDER_HEIGHT) )
    {
        printf( "Failed to initialize!\n" );
    }
    else
    {
        //try to load the texture
        textureParticle = loadTexture("../../assets/dot9.png", gRenderer);
        textureSphere = loadTexture("../../assets/sphere51.png", gRenderer);
        SDL_SetTextureColorMod(textureSphere, 128, 128, 128);
        if (textureParticle == NULL || textureSphere==NULL)
        {
            printf("Texture was not loaded. Exiting!\n");
            return 0;
        }

        //Main loop flag
        bool quit = false;
        bool pause = true;

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
                else if (e.type == SDL_KEYDOWN)
                {
                    if (e.key.keysym.sym == SDLK_p)
                     pause = !pause;
                }
            }

            // render_particles(*g_other, gRenderer);
            // SDL_Delay( 1000 ); //delay in millisecondes

            if (!pause)
            {
                std::cout << "iteration " << iter << std::endl;
                solver.solve(dt, iter);
            }
            
            SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
            SDL_RenderClear( gRenderer );

            render_spheres(solver, gRenderer, textureSphere);
            render_particles(solver, gRenderer, textureParticle);

            // SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
            // for (int i = 0; i < 20; i++)
            //     SDL_RenderDrawPoint( gRenderer, x+i+50, y+i+50 );
            // SDL_Rect srcrect = {0,0,9,9};
            // SDL_Rect dstrect = {0,0,20,20};
            // SDL_RenderCopy(gRenderer, gTexture, &srcrect, &dstrect);

            SDL_RenderPresent(gRenderer);

            SDL_Delay( 10 ); //delay in millisecondes

            //save the surface as png
            if (!pause){
                sprintf(filename_buf, "./pngs/frame_%04d.png", iter);

                SDL_RenderReadPixels(gRenderer, NULL, SDL_PIXELFORMAT_ARGB8888, frm->pixels, frm->pitch);
                IMG_SavePNG(frm, filename_buf);
                
                iter++;
            }
        }   
    }

    //Free resources and close SDL
    close(gWindow, gRenderer);
    SDL_FreeSurface(frm); 

    return 0;
}