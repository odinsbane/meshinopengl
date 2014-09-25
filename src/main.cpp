//
//  main.cpp
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//
#define GLFW_DISPLAY
#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#ifdef GLFW_DISPLAY
  #include "Display.h"
#endif
#include "rod.h"
#include "FRandom.h"
const double DT = 0.01;
const double rebound = 1;
std::vector<Rod*> rods;
void step();
void initialize();
void display();

#ifdef GLFW_DISPLAY
void initializeGraphics();
int updateGraphics();
Display* graphics;



#endif
int main(int argc, const char * argv[])
{
    std::cout<< "initializing\n";
    initialize();
#ifdef GLFW_DISPLAY
    initializeGraphics();
#endif
    
    bool run = true;
    int steps = 0;
	
	std::cout<< "running\n";

    while(run){
        
        for(int i = 0; i<5; i++){
            step();
        }
        steps++;
#ifdef GLFW_DISPLAY
        if(updateGraphics()!=0) run = false;
#else
        display();
        if(steps==10) run=false;
#endif
        
    }

    #ifdef GLFW_DISPLAY
        graphics->shutdown();
    #endif
    return 0;
}

#ifdef GLFW_DISPLAY
void initializeGraphics(){
    
    graphics = new Display(rods.size());
    graphics->initialize();

}

int updateGraphics(){
    int i = 0;
    for(auto itr = rods.begin(), end = rods.end(); itr!=end; itr++){
        Rod* b1 = *itr;
        graphics->updateRod(i, *b1);
        i++;
    }
    return graphics->render();
}
#endif



FRandom* generator;
void initialize(){
    generator = new FRandom(1);
    double delta = 0.2;
    for (int i = 0; i<20; i++){
        for(int j = 0; j<20; j++) {
            Rod *rod = new Rod(1.6, 0.075);
            rod->position[0] = delta*i - 1;
            rod->position[1] = delta*j - 1;
            glm::dvec3 dir = generator->randomDirection();

            rod->direction[0] = dir[0];
            rod->direction[1] = dir[1];
            rod->direction[2] = dir[2];
            rod->updateBounds();
            rods.push_back(rod);
        }
    }

}


void step(){
    int N = rods.size();
    for(int i = 0; i<N; i++){
        Rod *rod = rods[i];
        for(int j = i+1; j<N; j++){
            Rod *other = rods[j];
            double d = rod->collide(*other);
        }
    }

    for(int i = 0; i<N; i++){
        Rod *rod = rods[i];
        rod->prepareForces();
        rod->update(DT);
        rod->updateBounds();
        rod->clearForces();

    }
}

void display(){
    Rod* rod = rods[0];
    std::cout <<rod->position[0] << "\t" << rod->position[1] << "\n";
}
