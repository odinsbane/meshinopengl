
#ifndef ___PARALLELRODS__CONSTANTS__
#define ___PARALLELRODS__CONSTANTS__

#ifndef PI
#define PI 3.141592653589793
#endif

namespace Constants {

    //initialization
    const int ACTINS = 800;
    const int MYOSINS = 20;

    //simulation
    const double DT = 1e-3;
    const double WIDTH = 6.4;
    const double SEED_WIDTH=6.4;
    const double THICKNESS = 0.2;
    const double STEPS_PER_SIMULATE=1000;
    const double STEPS_PER_FRAME=1e2;
    const double SUB_STEPS=1;
    const double RELAXATION_LIMIT = 5e-1;
    const double ERROR_THRESHOLD = 1e-5;
    const double REPULSION=100;  //spring type force.
    const double MEMBRANE_POSITION=100;
    const double MEMBRANE_REPULSION = 0.1; //constant force
    const double BOTTOM_BOUNDS = 0.56;
    //parameters
    const double MYOSIN_LENGTH = 0.8; //300nm
    const double MYOSIN_DIAMETER = 0.2; //50nm
    const double MYOSIN_ACTIVE_FORCE = 1;
    const double MYOSIN_ALPHA_S = 1;
    const double MYOSIN_ALPHA = 1;
    const double K_m = 100.0;
    const double MYOSIN_BIND_LENGTH = 0.3;
    const double MYOSIN_BINDING_TIME = 10e32;
    const double MYOSIN_BIND_PROBABILITY = 0.9;
    const double MYOSIN_DIFFUSION_FORCE = 0;

    const double ACTIN_LENGTH = 2.0;  //500 nm.
    const double ACTIN_DIAMETER = 0.032;//8 nm.
    const double ACTIN_ALPHA = 1.0;
    const double ACTIN_LENGTH_SIGMA = 0.1;
    const double MEMBRANE_BOUND_ACTIN = 0.0;

    const double ANGLE_SIGMA = 0.0;
    const double CROSS_LINK_LENGTH = 0.2;
    const double CROSS_LINK_BINDING_TIME = 20;
    const double CROSS_LINK_BIND_PROBABILITY = 5e-3;
    const bool CROSS_LINKER_TURNOVER = false;
    const double K_x = 100;
    const bool STERIC_INTERACTIONS = false;

}

#endif
