

#ifndef __ParallelBalls__Simulation__
#define __ParallelBalls__Simulation__

#include <stdio.h>
#include <vector>
#include "rod.h"
#include "FRandom.h"
#include <math.h>
#include "interactions.h"
#include "Constants.h"
class Simulation{
    double working_dt;
    std::vector<std::unique_ptr<std::vector<double>>> position_record;
    std::vector<std::unique_ptr<std::vector<double>>> force_record;
    std::vector<ActinFilament*> actins;
    std::vector<MyosinMotor*> myosins;
    std::vector<MyosinMotorBinding*> bindings;
    std::vector<CrosslinkedFilaments*> xlinkers;
    FRandom* number_generator;
    void seedActinFilaments();
    void seedMyosinMotors();
    void freeSeedActinFilaments();
    void restorePositions(int index);
    void clearForces();
    void copyForces(int index);
    void copyPositions(int index);
    ActinFilament* createNewFilament();
    MyosinMotor* createNewMotor();
    CrosslinkedFilaments *createNewCrosslinkedFilaments();

    void prepareForUpdate(int con_cout, const std::vector<double> &coefficients);
    void prepareRelaxSpace();
    double prepareForces();
    void relax();
    void partialUpdate(double dt);
    void updateInteractions(double dt);
    double calculateError();
    glm::dvec3 getReflectedPoint(glm::dvec3 a, glm::dvec3 b);
    void applyMembraneForce(Rod* fil);
    public:
        Simulation(){ number_generator = new FRandom();}
        void initialize();
        void step();
        std::vector<ActinFilament*> &getActins();
        std::vector<MyosinMotor*> &getMyosins();

        void seedCrosslinkers();
        double getReflectedApproach(ActinFilament *a, ActinFilament *b);
        double reflectedCollision(Rod *other, Rod *filament);
        glm::dvec2 getReflectedIntersections(ActinFilament *a, ActinFilament *b);
        void crosslinkFilaments(ActinFilament *a, ActinFilament *b);
        void createTestCase();
        void twoFilamentTestCase();
        void seedCrosslinkerTestCase();
};

#endif
