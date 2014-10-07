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
#include <rod.h>

class MyosinMotorBinding{
    MyosinMotor& motor;
    std::array<ActinFilament*, 2> bound;

    public:
        static int FRONT=0;
        static int BACK = 0;

        MyosinMotorBinding(MyosinMotor &m){
            motor = m;
        }

};

class CrosslinkedFilaments{
    std::array<ActinFilament*, 2> filaments;
    std::array<double, 2> locations;

};

#endif