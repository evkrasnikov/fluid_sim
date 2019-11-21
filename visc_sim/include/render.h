#pragma once

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <string>

//Screen dimension constants
//const int SCREEN_WIDTH = 640;
//const int SCREEN_HEIGHT = 480;

//Starts up SDL and creates window
bool init(SDL_Window* &gWindow, SDL_Renderer* &gRenderer, int screen_width, int screen_height);

//Loads media
bool loadMedia(SDL_Surface* &gHelloWorld);

SDL_Texture* loadTexture( std::string path, SDL_Renderer* &gRenderer );

//Frees media and shuts down SDL
void close(SDL_Window* &gWindow, SDL_Renderer* &gRenderer);

