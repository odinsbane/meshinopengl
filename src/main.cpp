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
#include <chrono>
#include <QtCore/qlist.h>
#include <QtGui/qpolygon.h>

#ifdef GLFW_DISPLAY
  #include "Display.h"
#endif

#include "rod.h"

std::vector<Rod*> rods;
Simulation sim;

#ifdef GLFW_DISPLAY
    void initializeGraphics();
    int updateGraphics();
    Display* graphics;
#endif


int main(int argc, const char * argv[])
{
    std::cout<< "initializing\n";

    sim.initialize();

#ifdef GLFW_DISPLAY
    std::vector<ActinFilament*> &actins = sim.getActins();
    for(int i = 0;i<Constants::ACTINS; i++){
        Rod* a = actins[i];
        rods.push_back(a);
    }

    std::vector<MyosinMotor*> &myosins = sim.getMyosins();

    for(int i = 0; i<Constants::MYOSINS; i++){
        Rod* a = myosins[i];
        rods.push_back(a);
    }

    printf("Rods loaded: %ld\n", rods.size());
    initializeGraphics();
    updateGraphics();
#endif
    std::thread main([](){
        bool run = true;
        int steps = 0;

        std::cout<< "running\n";

        while(run){
            auto start = std::chrono::system_clock::now();
            for(int i = 0; i<Constants::STEPS_PER_FRAME; i++){

                sim.step();

            }
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            printf("step: %lld \n", elapsed);

            steps++;
            run = steps<=Constants::STEPS_PER_SIMULATE;

    #ifdef GLFW_DISPLAY
            if(updateGraphics()!=0) run = false;
    #else
            display();
            if(steps==10) run=false;
    #endif

        }

    });
    graphics->graphicsLoop();
    main.join();
    #ifdef GLFW_DISPLAY
        graphics->shutdown();
    #endif
    return 0;
}

#ifdef GLFW_DISPLAY
void initializeGraphics(){
    
    graphics = new Display(rods.size());
    graphics->initialize();
    //graphics->startGraphicsLoop();
}

int updateGraphics(){
    //graphics->render();
    int i = 0;
    for(auto itr = rods.begin(), end = rods.end(); itr!=end; itr++){
        Rod* b1 = *itr;
        graphics->updateRod(i, *b1);
        i++;
    }

    int springs = 0;
    std::vector<CrosslinkedFilaments*> &linked = sim.getCrosslinkedFilaments();
    std::vector<MyosinMotorBinding*> &bindings = sim.getMyosinMotorBindings();
    graphics->setSpringCount(linked.size() + 2*bindings.size());
    int index = 0;
    glm::dvec3 a,b,c;
    for(CrosslinkedFilaments* x: linked){
        a = x->filaments[0]->getPoint(x->locations[0]);
        b = x->filaments[1]->getPoint(x->locations[1]);
        c = sim.getReflectedPoint(a, b);
        graphics->updateSpring(index, a, c);
        index++;

    }
    int h;
    for(MyosinMotorBinding* bind: bindings){
        h = MyosinMotor::FRONT;

        a = bind->getHeadPosition(h);
        b = bind->getBindingPosition(h);
        c = sim.getReflectedPoint(a,b);
        graphics->updateSpring(index, a, c);
        index++;

        h = MyosinMotor::BACK;

        a = bind->getHeadPosition(h);
        b = bind->getBindingPosition(h);
        c = sim.getReflectedPoint(a,b);
        graphics->updateSpring(index, a, c);
        index++;



    }



    return graphics->getRunning();
}
#endif