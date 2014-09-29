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
const double DT = 0.1;
const double rebound = 1;
std::vector<Rod*> rods;
void step();
void display();

#ifdef GLFW_DISPLAY
    void initializeGraphics();
    int updateGraphics();
    Display* graphics;
#endif

int main(int argc, const char * argv[])
{
    std::cout<< "initializing\n";
    Simulation sim;
    sim.initialize();

#ifdef GLFW_DISPLAY
    std::vector<Rod*> &actins = sim.getActins();
    printf("%d\n", actins.size());
    for(int i = 0;i<Constants::ACTINS; i++){
        Rod* a = actins[i];

        rods.push_back(actins[i]);

    }
    printf("%d\n", rods.size());
    initializeGraphics();
#endif
    
    bool run = true;
    int steps = 0;
	
	std::cout<< "running\n";

    while(run){
        
        for(int i = 0; i<5; i++){

            sim.step();

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


void display(){
    Rod* rod = rods[0];
    std::cout <<rod->position[0] << "\t" << rod->position[1] << "\n";
}
