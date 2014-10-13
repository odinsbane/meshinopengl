
#ifndef ___PARALLELRODS__CONSTANTS__
#define ___PARALLELRODS__CONSTANTS__
namespace Constants {

    //initialization
    const int ACTINS = 400;
    const int MYOSINS = 100;

    //simulation
    const double DT = 1e-2;
    const double WIDTH = 6.4;
    const double SEED_WIDTH=6.4;
    const double THICKNESS = 0.56;
    const double STEPS_PER_SIMULATE=250;
    const double STEPS_PER_FRAME=1000;
    const double SUB_STEPS=10000;
    const double RELAXATION_LIMIT = 1;
    const double ERROR_THRESHOLD = 1e-5;
    const double REPULSION=0;  //spring type force.
    const double MEMBRANE_POSITION=0.4;
    const double MEMBRANE_REPULSION = 0.1; //constant force
    const double BOTTOM_BOUNDS = 0.56;
    //parameters
    const double MYOSIN_LENGTH = 0.8; //300nm
    const double MYOSIN_DIAMETER = 0.2; //50nm
    const double MYOSIN_ACTIVE_FORCE = 1;
    const double MYOSIN_ALPHA_S = 1;
    const double MYOSIN_ALPHA = 1;
    const double K_m = 1.0;
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

#endif