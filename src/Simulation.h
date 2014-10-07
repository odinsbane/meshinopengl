

#ifndef __ParallelBalls__Simulation__
#define __ParallelBalls__Simulation__

#include <stdio.h>
#include <vector>
#include "rod.h"
#include "FRandom.h"
#include <math.h>

namespace Constants {

    //initialization
    const int ACTINS = 240;
    const int MYOSINS = 80;

    //simulation
    const double DT = 1e-2;
    const double WIDTH = 6.4;
    const double SEED_WIDTH=6.4;
    const double THICKNESS = 0.56;
    const double STEPS_PER_SIMULATE=250;
    const double STEPS_PER_FRAME=1000;
    const double SUB_STEPS=1;
    const double RELAXATION_LIMIT = 500;
    const double ERROR_THRESHOLD = 1e+1;
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
    double working_dt;
    std::vector<std::unique_ptr<std::vector<double>>> position_record;
    std::vector<std::unique_ptr<std::vector<double>>> force_record;
    std::vector<ActinFilament*> actins;
    std::vector<MyosinMotor*> myosins;
    FRandom* number_generator;
    void seedActinFilaments();
    void seedMyosinMotors();
    void restorePositions(int index);
    void clearForces();
    void copyForces(int index);
    void copyPositions(int index);
    ActinFilament* createNewFilament();
    MyosinMotor* createNewMotor();
    void prepareForUpdate(int con_cout, const std::vector<double> &coefficients);
    void prepareRelaxSpace();
    double prepareForces();
    void relax();
    void partialUpdate(double dt);
    double calculateError();
    glm::dvec3 getReflectedPoint(glm::dvec3 a, glm::dvec3 b);
    public:
        Simulation(){ number_generator = new FRandom();}
        void initialize();
        void step();
        std::vector<ActinFilament*> &getActins();
        std::vector<MyosinMotor*> &getMyosins();

};

#endif
