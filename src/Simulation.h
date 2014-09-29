

#ifndef __ParallelBalls__Simulation__
#define __ParallelBalls__Simulation__

#include <stdio.h>
#include <vector>
#include "rod.h"
#include "FRandom.h"
#include <math.h>

namespace Constants {

    //initialization
    const int ACTINS = 1000;
    const int MYOSINS = 0;

    //simulation
    const double DT = 1e-1;
    const double WIDTH = 6.4;
    const double SEED_WIDTH=6.4;
    const double THICKNESS = 0.56;
    const double STEPS_PER_SIMULATE=250;
    const double STEPS_PER_FRAME=1000;
    const double SUB_STEPS=10000;
    const double RELAXATION_LIMIT = 5e-1;
    const double ERROR_THRESHOLD = 1e-12;
    const double REPULSION=1;  //spring type force.
    const double MEMBRANE_POSITION=0.4;
    const double MEMBRANE_REPULSION = 1; //constant force
    const double BOTTOM_BOUNDS = 0.56;
    //parameters
    const double MYOSIN_LENGTH = 0.8; //300nm
    const double MYOSIN_DIAMETER = 0.2; //50nm
    const double MYOSIN_ACTIVE_FORCE = 1;
    const double MYOSIN_ALPHA_S = 1;
    const double MYOSIN_ALPHA = 1;
    const double K_m = 100;
    const double MYOSIN_BIND_LENGTH = 0.2;
    const double MYOSIN_BINDING_TIME = 10e32;
    const double MYOSIN_BIND_PROBABILITY = 0.1;
    const double MYOSIN_DIFFUSION_FORCE = 0;

    const double ACTIN_LENGTH = 2.0;  //500 nm.
    const double ACTIN_DIAMETER = 0.032;//8 nm.
    const double ACTIN_ALPHA = 0.1;
    const double ACTIN_LENGTH_SIGMA = 0.1;
    const double MEMBRANE_BOUND_ACTIN = 0.0;

    const double CROSS_LINK_LENGTH = 0.2;
    const double CROSS_LINK_BINDING_TIME = 10000;
    const double CROSS_LINK_BIND_PROBABILITY = 0.1;
    const bool CROSS_LINKER_TURNOVER = false;
    const double K_x = 100;

}

class Simulation{
    std::vector<double*> position_record;
    std::vector<double*> torque_record;
    std::vector<double*> force_record;
    std::vector<Rod*> actins;
    std::vector<Rod*> myosins;
    FRandom* number_generator;
    void seedActinFilaments();
    void seedMyosinMotors();
    public:
        Simulation(){ number_generator = new FRandom();}
        void initialize();
        void step();
        std::vector<Rod*> &getActins();


};

#endif
