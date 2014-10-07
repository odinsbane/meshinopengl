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
class MyosinMotorBinding{
    MyosinMotor* motor;
    std::array<double, 2> binding_position;
    std::array<double, 2> current_time;
    std::array<double, 2> unbind_time;
    FRandom* number_generator;

    public:


        MyosinMotorBinding(MyosinMotor* m){
            motor = m;
        }
        void bind(ActinFilament* f, int head, double position);
        void applyForces();
        void setNumberGenerator(FRandom* ng){number_generator=ng;}
};

class CrosslinkedFilaments{
    std::array<ActinFilament*, 2> filaments;
    std::array<double, 2> locations;

};

#endif