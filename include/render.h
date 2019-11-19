#pragma once

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

//Screen dimension constants
const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;

//Starts up SDL and creates window
bool init(SDL_Window* &gWindow, SDL_Renderer* &gRenderer);

//Loads media
bool loadMedia(SDL_Surface* &gHelloWorld);

//Frees media and shuts down SDL
void close(SDL_Window* &gWindow, SDL_Renderer* &gRenderer);

