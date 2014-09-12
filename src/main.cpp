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
#include "Ball.h"
#include <math.h>
#ifdef GLFW_DISPLAY
  #include "Display.h"
#endif
#include <functional>
#include "ExecutionService.h"

const double DT = 0.01;
const double rebound = 1;
std::vector<Ball*> balls;
ExecutionService service(32);

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
    
    graphics = new Display(balls.size());
    graphics->initialize();

}

int updateGraphics(){
    int i = 0;
    for(auto itr = balls.begin(), end = balls.end(); itr!=end; itr++){
        Ball* b1 = *itr;
        graphics->updateBall(i, b1->X, b1->Y, b1->radius);
        i++;
    }
    return graphics->render();
}
#endif

void initialize(){
    
    double delta = 0.014;
    for (int i = 0; i<75; i++){
        for(int j = 0; j<140; j++){
            Ball* ball = new Ball(i*delta - 0.2, j*delta - 0.97);
            balls.push_back(ball);
            
        }
    }
    Ball* ball = new Ball(-0.5, 0);
    ball->Vx = 0.5;
    ball->Vy = 0.01;
    balls.push_back(ball);


}


void step(){
    for(auto itr = balls.begin(), end = balls.end(); itr!=end; itr++){
        auto f = [itr, end](){
            Ball* b1 = *itr;
            for(auto otra = itr+1; otra!=end; otra++){
                Ball* b2 =*otra;
                double min = b1->radius + b2->radius;
                min = min*min;
                double dx = b1->X - b2->X;
                if(dx*dx>min) continue;
                
                double dy = b1->Y - b2->Y;
                if(dy*dy>min) continue;
                
                double mag = dx*dx + dy*dy;
                if(mag<min){
                    mag = sqrt(mag);
                    dx = dx/mag;
                    dy = dy/mag;

                    double k = (b2->Vx - b1->Vx)*dx + (b2->Vy - b1->Vy)*dy;
                    double fx = k*dx/DT;
                    double fy = k*dy/DT;
                    double dot = fx*dx + fy*dy;
                    if(dot>0) {
                        b2->applyForce(-fx, -fy);
                        b1->applyForce(fx, fy);
                    }
                }
            }
            
            
            //b1->applyForce(0, -0.01);
        };
        service.submit(f);
    }
    
    service.waitForExecution();
    
    for(auto itr = balls.begin(), end = balls.end(); itr!=end; itr++){
        Ball* b1 = *itr;
        service.submit([b1](){
            b1->update(DT);
            if(b1->X - b1->radius < -1.0){
                b1->X = b1->radius - 1;
                b1->Vx = fabs(b1->Vx);
            } else if(b1->X + b1->radius > 1.0){
                b1->X = -b1->radius + 1;
                b1->Vx = -fabs(b1->Vx);
            }
            
            if(b1->Y - b1->radius<-1){
                b1->Y = b1->radius - 1;
                b1->Vy = fabs(b1->Vy);
            } else if(b1->Y + b1->radius>1){
                b1->Y = 1 - b1->radius;
                b1->Vy = -fabs(b1->Vy);
            }
        });
    }
    
    service.waitForExecution();

}

void display(){
    Ball* ball = balls[0];
    std::cout <<ball->X << "\t" << ball->Y << "\n";
}
