//
//  Ball.cpp
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#include "Ball.h"

Ball::Ball(double x, double y){
    X = x;
    Y = y;
    Vx = 0;
    Vy = 0;
    radius = 0.006;
}

void Ball::applyForce(double fx, double fy){
    double* force = new double[2];
    force[0] = fx;
    force[1] = fy;
    //std::shared_ptr<double> fptr(force, [](double* d){delete[] d;});
    std::lock_guard<std::mutex> lock(mutex);
    forces.push_back(force);
}

void Ball::update(double dt){
    
    std::lock_guard<std::mutex> lock(mutex);
    
    double fx = 0;
    double fy = 0;
    
    for(auto itr = forces.begin(), end = forces.end(); itr!=end; itr++){
        double* force = *itr;
        fx += force[0];
        fy += force[1];
        delete force;
    }
    
    Vx += fx*dt - 0.02*Vx*dt;
    Vy += fy*dt - 0.02*Vy*dt;
    
    X += Vx*dt;
    Y += Vy*dt;
    
    forces.clear();
}
