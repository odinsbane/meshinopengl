#include "Simulation.h"

Rod* createNewFilament(){
    Rod* a = new Rod(Constants::ACTIN_LENGTH, Constants::ACTIN_DIAMETER);
    a->alpha_longitudinal = Constants::ACTIN_ALPHA;
    a->alpha_perpendicular = Constants::ACTIN_ALPHA;
    a->alpha_rotational = Constants::ACTIN_ALPHA;
    return a;
}

const double PI = 3.141592653589793;


void Simulation::seedActinFilaments(){
    for(int i = 0; i<Constants::ACTINS; i++){
        if(number_generator->nextDouble()<Constants::MEMBRANE_BOUND_ACTIN) {
            //generate a bound filament.
            double x = Constants::SEED_WIDTH * number_generator->nextDouble() - 0.5 * Constants::SEED_WIDTH;
            double y = Constants::SEED_WIDTH * number_generator->nextDouble() - 0.5 * Constants::SEED_WIDTH;

            Rod* f = createNewFilament();
            //Crosslinker link = createNewCrossLinker();

            double min = Constants::MEMBRANE_POSITION - 0.5*Constants::THICKNESS;
            double max = Constants::MEMBRANE_POSITION + 0.5*Constants::THICKNESS;

            double l = f->length*0.5;// + link.length; //from center of filament to xlink point.
            max = l>max?max:l;
            if(min>l){
                printf("Filaments are too short to bind membrane\n");
                exit(-1);
            }
            double z = number_generator->nextDouble()*(max - min) + min;
            double r = sqrt(l*l - z*z);

            double theta = PI*2*number_generator->nextDouble();
            double dx = cos(theta)*r;
            double dy = sin(theta)*r;

            f->position = glm::dvec3(
                    x - dx,
                    y - dy,
                    Constants::MEMBRANE_POSITION - z
            );
            glm::dvec3 pt(x, y, Constants::MEMBRANE_POSITION);
            f->direction = glm::normalize(pt - f->position);

            // SlidingMembraneCrossLinker linkage = new SlidingMembraneCrossLinker(this, f, link, 1e32);

            //addCrossLinking(linkage);
            //addXLinker(link);
            actins.push_back(f);
            //relaxRod(f);
        } else {
            //generate a free filament.
            double x = Constants::SEED_WIDTH * number_generator->nextDouble() - 0.5 * Constants::SEED_WIDTH;
            double y = Constants::SEED_WIDTH * number_generator->nextDouble() - 0.5 * Constants::SEED_WIDTH;
            double z = Constants::THICKNESS * number_generator->nextDouble() - 0.5 * Constants::THICKNESS;

            double theta = 2 * PI * number_generator->nextDouble();


            double distance_to_membrane = Constants::MEMBRANE_POSITION - z;
            Rod* f = createNewFilament();

            double available_angle = distance_to_membrane > (0.5 * f->length) ? PI : asin(distance_to_membrane / (0.5 * f->length));
            double phi = PI / 2 + (0.5 - number_generator->nextDouble()) * available_angle;
            f->direction = glm::dvec3(
                    cos(theta) * sin(phi),
                    sin(theta) * sin(phi),
                    cos(phi)
            );


            f->position = glm::dvec3(
                    x,
                    y,
                    z
            );
            f->updateBounds();
            actins.push_back(f);
            //relaxRod(f);
        }
    }
}

void Simulation::prepareRelaxSpace(){
    working_dt = Constants::DT;
    std::vector<double> *ps = new std::vector<double>(6*actins.size() + 6*myosins.size());
    position_record.push_back(std::unique_ptr<std::vector<double>>(ps));

}
void Simulation::seedMyosinMotors(){

}

void Simulation::initialize(){

    seedActinFilaments();
    seedMyosinMotors();

}

void Simulation::prepareForces(){
    int N = actins.size();
    for(int i = 0; i<N; i++){
        Rod *rod = actins[i];
        for(int j = i+1; j<N; j++){
            Rod *other = actins[j];
            double d = rod->collide(*other);
        }
    }

    for(int i = 0; i<N; i++){
        Rod *rod = actins[i];
        rod->prepareForces();
    }


}

void Simulation::copyPositions(int index){
    int i = 0;
    std::vector<double> *positions = position_record[index].get();
    for(Rod* & rod : actins){
        for(int j = 0; j<3; j++){
            positions[0][i + j]= rod->position[j];
            positions[0][i + j + 3] = rod->direction[j];
        }
        i+=6;
    }

    for(Rod* & rod : myosins){
        for(int j = 0; j<3; j++){
            positions[0][i + j] = rod->position[j];
            positions[0][i + j + 3] = rod->direction[j];
        }
        i+=6;
    }

}

void Simulation::restorePositions(int index){
    int i = 0;
    auto positions = position_record[index].get();
    for(Rod* & rod : actins){
        for(int j = 0; j<3; j++){
            (*positions)[i + j]= rod->position[j];
            (*positions)[i + j + 3] = rod->direction[j];
        }
        i+=6;
    }

    for(Rod* & rod : myosins){
        for(int j = 0; j<3; j++){
            (*positions)[i + j] = rod->position[j];
            (*positions)[i + j + 3] = rod->direction[j];
        }
        i+=6;
    }
}

void Simulation::copyForces(int index){
    int i = 0;
    std::vector<double> *forces = force_record[index].get();
    for(Rod* & rod : actins){
        for(int j = 0; j<3; j++){
            forces[0][i + j]= rod->force[j];
            forces[0][i + j + 3] = rod->torque[j];
        }
        i+=6;
    }

    for(Rod* & rod : myosins){
        for(int j = 0; j<3; j++){
            forces[0][i + j] = rod->force[j];
            forces[0][i + j + 3] = rod->torque[j];
        }
        i+=6;
    }
}

void Simulation::clearForces(){

    for(Rod* & rod : actins){
        rod->clearForces();
    }

    for(Rod* & rod : myosins){
        rod->clearForces();
    }
}

void Simulation::prepareForUpdate(int con_count, const std::vector<double> &coefficients) {
    int i = 0;
    for(Rod* & rod : actins){
        for(int k = 0; k<con_count; k++){
            std::vector<double> *forces = force_record[k].get();
            for(int j = 0; j<3; j++) {
                rod->force[j] += coefficients[k] * forces[0][i + j];
                rod->torque[j] += coefficients[k] * forces[0][i + j + 3];
            }
        }
        i+=6;
    }

    for(Rod* & rod : myosins){
        for(int k = 0; k<con_count; k++){
            std::vector<double> *forces = force_record[k].get();
            for(int j = 0; j<3; j++) {
                rod->force[j] += coefficients[k] * forces[0][i + j];
                rod->torque[j] += coefficients[k] * forces[0][i + j + 3];
            }
        }

        i+=6;
    }
}

void Simulation::partialUpdate(double dt){
    for(Rod* & rod : actins){
        rod->update(dt);
        rod->clearForces();
    }

    for(Rod* & rod : myosins){
        rod->update(dt);
        rod->clearForces();
    }
}
/**
* Using the Runge-Kutta-Fehlberg Method RKF45 an adaptive timestep technique
*
*/
void Simulation::relax(){
    //copy yk into position and direction record
    copyPositions(0);

    //create k1/h
    prepareForces();
    copyForces(0);
    clearForces();

    //set the forces to be 1/(4h) k1
    double args[1]{0.25};
    prepareForUpdate(1, {0.25});

    //update the position to y0 + 1/4k1
    partialUpdate(working_dt);

    //prepare k2
    prepareForces();
    copyForces(1);
    clearForces();
    //move back to y0
    restorePositions(0);
    prepareForUpdate(2, {3.0/32.0, 9.0/32.0});

    //y = y0 + 3/32k1 + 9/32k2
    partialUpdate(working_dt);

    //k3
    prepareForces();
    copyForces(2);
    clearForces();
    restorePositions(0);

    prepareForUpdate(3, {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0});
    //y = y0 + 1932/2197 k1 - 7200/2197 k2 + 7296/2197 k3
    partialUpdate(working_dt);

    //k4
    prepareForces();
    copyForces(3);
    clearForces();
    restorePositions(0);

    prepareForUpdate(4, {439.0/216.0, -8, 3680.0/513.0, -845.0/4104.0});
    //y = y0 + 439/216 k1 - 8 k2 + 3680/513 k3 - 845/4104 k4
    partialUpdate(working_dt);

    //k5
    prepareForces();
    copyForces(4);
    clearForces();
    restorePositions(0);

    prepareForUpdate(5, {-8.0/27.0, 2, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0});
    //y =y0 -8/27 k1 + 2 k2 - 3544/2565 k3 + 1859/4104 k4 - 11/40 k5
    partialUpdate(working_dt);


}



void Simulation::step(){
    int N = actins.size();
    for(int i = 0; i<N; i++){
        Rod *rod = actins[i];
        for(int j = i+1; j<N; j++){
            Rod *other = actins[j];
            double d = rod->collide(*other);
        }
    }

    for(int i = 0; i<N; i++){
        Rod *rod = actins[i];
        rod->prepareForces();
        rod->update(Constants::DT);
        rod->updateBounds();
        rod->clearForces();

    }
}

std::vector<Rod*> &Simulation::getActins(){
    return actins;
}