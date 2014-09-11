//
//  Ball.h
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#ifndef __ParallelBalls__Ball__
#define __ParallelBalls__Ball__

#include <iostream>
#include <vector>
#include <mutex>

class Ball{
private:
    std::vector<double*> forces;
    std::mutex mutex;
public:
    double X,Y, Vx, Vy, radius;
    void update(double dt);
    void applyForce(double fx, double fy);
    Ball(double x, double y);
    ~Ball(){}
};

#endif /* defined(__ParallelBalls__Ball__) */
