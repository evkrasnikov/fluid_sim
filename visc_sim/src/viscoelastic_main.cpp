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
    
    for(int k = 0; k < solver.particles.size(); k++)
    {
        Particle& p = solver.particles[k];
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
        Eigen::Vector2d particle_pos = p.pos;
        Eigen::Vector2d world_pos = s*(particle_pos-c);

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

    for(int k = 0; k < solver.spheres.size(); k++)
    {
        Sphere s = solver.spheres[k];

        Eigen::Vector2d sphere_pos = s.pos;
        Eigen::Vector2d world_pos = s_mat*(sphere_pos-c);
        Eigen::Vector2d world_r; 
        world_r << s.r*RENDER_WIDTH/WORLD_WIDTH, s.r*RENDER_HEIGHT/WORLD_HEIGHT;

        int x = round(world_pos(0)-world_r(0));
        int y = round(world_pos(1)-world_r(1));
        SDL_Rect srcrect = {0,0,201,201};
        SDL_Rect dstrect = {x,y,int(world_r(0)*2),int(world_r(1)*2)};
        SDL_RenderCopy(gRenderer, texture, &srcrect, &dstrect);
    }
}

// zero gravity,2 water blobs merging together
void scene1(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;
    set.RHO0 = 15;
    set.K = 0.1;
    set.KNEAR = 1;
    set.SIGMA = 0;
    set.BETA  = 0.2; 
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 30; 
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 5;
    set.G << 0, 0;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 0; 

    WaterRect water_blob1 = {30, 30, 10, 80, 1, 0, 0};
    WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    water_blobs.push_back(water_blob2);
}

// zero gravity 2 water blobs eventually taking a circular shape
void scene2(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;
    set.RHO0 = 15;
    set.K = 25;
    set.KNEAR = 250;
    set.SIGMA = 0;
    set.BETA  = 0.2;
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 30; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 5;
    set.G << 0, 0;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 0; 

    WaterRect water_blob1 = {30, 30, 10, 80, 1, 0, 0};
    WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    water_blobs.push_back(water_blob2);
}

// different resolution (low resolution / large distance between the particles)
void scene3(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 5;
    set.RHO0 = 2;
    set.K = 0.5;
    set.KNEAR = 5;
    set.SIGMA = 0;
    set.BETA  = 0.2; 
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 30;
    set.D_STICK = 1.55;
    set.MU = 0.9;
    set.K_STICK = 0.5;//35;
    set.G << 0, -9.81;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 0; 

    WaterRect water_blob1 = {30, 40, 10, 80, 2, 0, 0};
    WaterRect water_blob2 = {30, 40, 60, 50, 2, 0, 0};
    water_blobs.push_back(water_blob1);
    water_blobs.push_back(water_blob2);
}

// different resolution (low resolution / small distance between the particles)
void scene4(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;
    set.RHO0 = 20;
    set.K = 0.5;
    set.KNEAR = 5;
    set.SIGMA = 0;
    set.BETA  = 0.2; 
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 30; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 0.5;//35;
    set.G << 0, -9.81;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 0; 

    WaterRect water_blob1 = {30, 40, 10, 80, 1, 0, 0};
    WaterRect water_blob2 = {30, 40, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    water_blobs.push_back(water_blob2);
}

// viscosity (low viscosity)
void scene5(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;//3;//
    set.RHO0 = 15;//10;//
    set.K = 0.5;//0.04;//
    set.KNEAR = 5;//0.1; //
    set.SIGMA = 0;//1;//
    set.BETA  = 0.03; //0.2;
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 30; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 5;//35;
    set.G << 0, -9.81;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 0; 

    WaterRect water_blob1 = {30, 30, 10, 100, 1, 0, 0};
    WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    water_blobs.push_back(water_blob2);
}

// viscosity (high viscosity)
void scene6(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;//3;//
    set.RHO0 = 15;//10;//
    set.K = 0.5;//0.04;//
    set.KNEAR = 5;//0.1; //
    set.SIGMA = 0.5;//1;//
    set.BETA  = 0.6; //0.2;
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 30; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 5;//35;
    set.G << 0, -9.81;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 0; 

    WaterRect water_blob1 = {30, 30, 10, 100, 1, 0, 0};
    WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    water_blobs.push_back(water_blob2);
}

// elastic and plastic material (high elasticity and plasticity)
void scene7(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;//3;//
    set.RHO0 = 15;//10;//
    set.K = 0.5;//0.04;//
    set.KNEAR = 5;//0.1; //
    set.SIGMA = 0;//1;//
    set.BETA  = 0.2; //0.2;
    set.GAMMA = 0.1;
    set.ALPHA = 0.5;
    set.K_SPRING = 20; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 25;//35;

    set.G << 0, -9.81;
    set.EN_MOLDING = 1;
    set.EN_SPRINGS = 1; 

    WaterRect water_blob1 = {80, 20, 10, 100, 1, 0, 0};
    //WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    //water_blobs.push_back(water_blob2);

    for (int i = 0; i < 4; i++)
    {
        Sphere s;
        s.pos << 10+i*20,50;
        s.r = 8;
        spheres.push_back(s);
    }
}

// elastic and plastic material (low elasticity and plasticity)
void scene8(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;//3;//
    set.RHO0 = 15;//10;//
    set.K = 0.5;//0.04;//
    set.KNEAR = 5;//0.1; //
    set.SIGMA = 0;//1;//
    set.BETA  = 0.2; //0.2;
    set.GAMMA = 0.2;
    set.ALPHA = 0.1;
    set.K_SPRING = 40; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 25;//35;

    set.G << 0, -9.81;
    set.EN_MOLDING = 1;
    set.EN_SPRINGS = 1; 

    WaterRect water_blob1 = {80, 20, 10, 100, 1, 0, 0};
    //WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    //water_blobs.push_back(water_blob2);

    for (int i = 0; i < 4; i++)
    {
        Sphere s;
        s.pos << 10+i*20,50;
        s.r = 8;
        spheres.push_back(s);
    }
}

// zero gravity,2 water blobs merging together
void scene9(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;//3;//
    set.RHO0 = 15;//10;//
    set.K = 5;//0.04;//
    set.KNEAR = 50;//0.1; //
    set.SIGMA = 0;//1;//
    set.BETA  = 0.2; //0.2;
    set.GAMMA = 0.1;
    set.ALPHA = 0.3;
    set.K_SPRING = 30; //0.3;
    set.D_STICK = 1.55;
    set.MU = 0.8;
    set.K_STICK = 5;//35;
    set.G << 0, 0;//-9.81;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 0; 

    WaterRect water_blob1 = {30, 30, 30, 80, 1, 0, 0};
    WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    water_blobs.push_back(water_blob2);
}

// zero gravity,2 water blobs merging together
void scene10(Settings& set, std::vector<Sphere>& spheres, std::vector<WaterRect>& water_blobs)
{
    set.H = 3.4;//3;//
    set.RHO0 = 15;//10;//
    set.K = 1.5;//0.04;//
    set.KNEAR = 15;//0.1; //
    set.SIGMA = 0;//1;//
    set.BETA  = 0.2; //0.2;
    
    // plasticity and elasticity 
    set.GAMMA = 0.1;
    set.ALPHA = 0.5;
    set.K_SPRING = 0.3; //0.3;
    set.D_STICK = 2.55;
    set.MU = 0.9;

    // boundary stuff
    set.K_STICK = 45;//35;

    set.G << 0, -9.81;
    set.EN_MOLDING = 0;
    set.EN_SPRINGS = 1; 

    WaterRect water_blob1 = {50, 30, 30, 110, 1, 0, 0};
    //WaterRect water_blob2 = {30, 30, 60, 50, 1, 0, 0};
    water_blobs.push_back(water_blob1);
    //water_blobs.push_back(water_blob2);

    Sphere s;
    s.pos << 50,50;
    s.r = 20;
    spheres.push_back(s);
    s.pos << 20, 20;
    s.r = 10;
    spheres.push_back(s);
}



int main(int argc, char* args[])
{
    if (argc != 3)
    {
        printf("Usage: <executable> scene_number max_iterations\n");
        return 0;
    }
    printf("%d\n", argc);
    int scene_num = atoi(args[1]);
    int max_iters = atoi(args[2]);
    printf("%d %d\n", scene_num, max_iters);

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
    //SDL_Surface *frm = SDL_CreateRGBSurfaceWithFormat(0, RENDER_WIDTH, RENDER_HEIGHT, 24, SDL_PIXELFORMAT_RGB24);

    std::vector<WaterRect> water_blobs;
    std::vector<Sphere> spheres;
    Settings set;

    switch(scene_num)
    {
        case 1: scene1(set, spheres, water_blobs); break;
        case 2: scene2(set, spheres, water_blobs); break;
        case 3: scene3(set, spheres, water_blobs); break;
        case 4: scene4(set, spheres, water_blobs); break;
        case 5: scene5(set, spheres, water_blobs); break;
        case 6: scene6(set, spheres, water_blobs); break;
        case 7: scene7(set, spheres, water_blobs); break;
        case 8: scene8(set, spheres, water_blobs); break;
        case 9: scene9(set, spheres, water_blobs); break;
        case 10: scene10(set, spheres, water_blobs); break;
        default: scene10(set, spheres, water_blobs); break;
    }

    FluidSolver solver(WORLD_WIDTH, WORLD_HEIGHT, set, spheres, water_blobs);
    double dt = 1.0/30.0;///4.0;

    //Start up SDL and create window
    if( !init(gWindow, gRenderer, RENDER_WIDTH, RENDER_HEIGHT) )
    {
        printf( "Failed to initialize!\n" );
    }
    else
    {
        //try to load the texture
        textureParticle = loadTexture("../../assets/dot9_white.png", gRenderer);
        textureSphere = loadTexture("../../assets/sphere201.png", gRenderer);
        SDL_SetTextureColorMod(textureParticle, 0, 128, 128);
        SDL_SetTextureColorMod(textureSphere, 128, 128, 128);
        if (textureParticle == NULL || textureSphere==NULL)
        {
            printf("Texture was not loaded. Exiting!\n");
            return 0;
        }

        //Main loop flag
        bool quit = false;
        bool pause = false;

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

            if (!pause)
            {
                std::cout << "iteration " << iter << std::endl;
                solver.solve(dt, iter);
            }
            
            SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
            SDL_RenderClear( gRenderer );

            render_spheres(solver, gRenderer, textureSphere);
            render_particles(solver, gRenderer, textureParticle);

            SDL_RenderPresent(gRenderer);

            SDL_Delay( 10 ); //delay in millisecondes

            //save the surface as png
            if (!pause){
                sprintf(filename_buf, "./pngs/frame_%04d.png", iter);

                SDL_RenderReadPixels(gRenderer, NULL, /*SDL_PIXELFORMAT_RGB24*/SDL_PIXELFORMAT_ARGB8888, frm->pixels, frm->pitch);
                IMG_SavePNG(frm, filename_buf);
                
                iter++;

                if (iter > max_iters)
                    break;
            }
        }   
    }

    //Free resources and close SDL
    close(gWindow, gRenderer);
    SDL_FreeSurface(frm); 

    return 0;
}