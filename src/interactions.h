//
//  interactions.h
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#ifndef __ParallelBalls__interactions__
#define __ParallelBalls__interactions__
#include <array>
#include "rod.h"
#include "FRandom.h"
#include <math.h>
#include <glm/geometric.hpp> //glm::dot
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include "Constants.h"

class MyosinMotorBinding{
    MyosinMotor* motor;
    std::array<double, 2> binding_position;
    std::array<double, 2> current_time;
    std::array<double, 2> unbind_time;
    std::array<double, 2> sliding;
    FRandom* number_generator;

    public:


        MyosinMotorBinding(MyosinMotor* m){
            motor = m;
        }
        void bind(ActinFilament* f, int head, double position);
        void applyForces();
        void setNumberGenerator(FRandom* ng){number_generator=ng;}
        void headForce(int head);

};

class CrosslinkedFilaments{
    std::array<double, 2> locations;
    double unbind_time;
    double current_time;
    public:
        std::array<ActinFilament*, 2> filaments;
        double K_x;
        double length;
        double tau_B;
        void applyForces();
        void update(double dt);
        void unbind();
};

#endif