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
#ifdef GLFW_DISPLAY
  #include "Display.h"
#endif
#include "rod.h"

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
        if(steps==100) run = false;
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

void initialize(){
    
    double delta = 0.014;
    for (int i = 0; i<21; i++){
        for(int j = 0; j<21; j++){
            Rod* rod = new Rod(0.1, 0.01);
            rod->direction[1] = 1;
            rod->position[0] = i*0.1 - 1;
            rod->position[1] = j*0.1 - 1;
            rods.push_back(rod);
            
        }
    }


}


void step(){

}

void display(){
    Rod* rod = rods[0];
    std::cout <<rod->position[0] << "\t" << rod->position[1] << "\n";
}
