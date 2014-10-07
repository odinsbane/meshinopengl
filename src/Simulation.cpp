#include "Simulation.h"

ActinFilament* Simulation::createNewFilament(){
    ActinFilament* a = new ActinFilament(Constants::ACTIN_LENGTH, Constants::ACTIN_DIAMETER);
    a->alpha_longitudinal = Constants::ACTIN_ALPHA;
    a->alpha_perpendicular = Constants::ACTIN_ALPHA;
    a->alpha_rotational = Constants::ACTIN_ALPHA;
    return a;
}

MyosinMotor* Simulation::createNewMotor(){
    MyosinMotor* m = new MyosinMotor(Constants::MYOSIN_LENGTH, Constants::MYOSIN_DIAMETER);
    m->F0 = Constants::MYOSIN_ACTIVE_FORCE;
    m->alpha_longitudinal = Constants::MYOSIN_ALPHA;
    m->alpha_perpendicular = Constants::MYOSIN_ALPHA;
    m->alpha_rotational = Constants::MYOSIN_ALPHA;
    m->alpha_s = Constants::MYOSIN_ALPHA_S;
    m->K_m = Constants::K_m;
    m->tau_B = Constants::MYOSIN_BINDING_TIME;

    return m;
}
const double PI = 3.141592653589793;


void Simulation::seedActinFilaments(){
    for(int i = 0; i<Constants::ACTINS; i++){
        if(number_generator->nextDouble()<Constants::MEMBRANE_BOUND_ACTIN) {
            //generate a bound filament.
            double x = Constants::SEED_WIDTH * number_generator->nextDouble() - 0.5 * Constants::SEED_WIDTH;
            double y = Constants::SEED_WIDTH * number_generator->nextDouble() - 0.5 * Constants::SEED_WIDTH;

            ActinFilament* f = createNewFilament();
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
            ActinFilament* f = createNewFilament();

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
    position_record.clear();
    force_record.clear();
    working_dt = Constants::DT;
    int n = 6 * Constants::ACTINS + 6 * Constants::MYOSINS;
    for(int i = 0; i<2; i++) {
        std::vector<double> *ps = new std::vector<double>(n);
        position_record.push_back(std::unique_ptr<std::vector<double>>(ps));
    }
    for(int i = 0; i<6; i++) {
        std::vector<double> *ps = new std::vector<double>(n);
        force_record.push_back(std::unique_ptr<std::vector<double>>(ps));
    }
}

glm::dvec3 Simulation::getReflectedPoint(glm::dvec3 a, glm::dvec3 b) {
    return b;
}
void Simulation::seedMyosinMotors(){
    for(int i = 0; i<Constants::MYOSINS; i++){
        MyosinMotor* motor = createNewMotor();
        std::cout << motor->isBound(0);
        MyosinMotorBinding* bind = new MyosinMotorBinding(motor);
        bind->setNumberGenerator(number_generator);
        ActinFilament* host = actins[number_generator->nextInt(actins.size())];

        double s = (number_generator->nextDouble() - 0.5) * host->length;
        glm::dvec3 host_a = host->getPoint(s);

        bind->bind(host, MyosinMotor::FRONT, s);
        std::cout << motor->isBound(0) <<"\n";
        std::vector<ActinFilament*> possibles;
        for (ActinFilament* target : actins) {
            if (host == target) continue;
            double separation = target->closestApproach(host_a);
            if (separation < motor->length) {
                glm::dvec3 reflected = getReflectedPoint(target->position, host_a);
                std::vector<double> intersections = target->getIntersections(reflected, motor->length + 2*Constants::MYOSIN_BIND_LENGTH);
                if (intersections.size() == 0) {
                    //The attached head of the motor is too close to the center of the target filament.
                    continue;
                }
                possibles.push_back(target);
            }
        }
        glm::dvec3 host_b;
        if (possibles.size() > 0) {

            ActinFilament* other = possibles[number_generator->nextInt(possibles.size())];
            glm::dvec3 reflected = getReflectedPoint(other->position, host_a);
            std::vector<double> intersections = other->getIntersections(reflected, motor->length + 2*Constants::MYOSIN_BIND_LENGTH);

            double s2 = intersections[number_generator->nextInt(intersections.size())];
            glm::dvec3 target_p = other->getPoint(s2);
            host_b = getReflectedPoint(host_a, target_p);
            bind->bind(other, MyosinMotor::BACK, s2);
        } else {
            do {
               double phi = PI * number_generator->nextDouble();
                double theta = 2 * PI * number_generator->nextDouble();
                //TODO remove
                double l = motor->length + motor->diameter + host->diameter;
                host_b = glm::dvec3(
                        host_a[0] + l * cos(theta) * sin(phi),
                        host_a[1] + l * sin(theta) * sin(phi),
                        host_a[2] + l * cos(phi)
                );
                printf("myosin head could not find a filament to attach to.\n");
            } while (host_b[2] > Constants::MEMBRANE_POSITION);

        }


        bindings.push_back(bind);


        motor->position = glm::dvec3(
                (host_a[0] + host_b[0]) * 0.5,
                (host_a[1] + host_b[1]) * 0.5,
                (host_a[2] + host_b[2]) * 0.5
        );


        glm::dvec3 dir = host_a - host_b;

        motor->direction = glm::normalize(dir);

        motor->updateBounds();

        myosins.push_back(motor);
        //relaxRod(motor);
    }
}

void Simulation::initialize(){
    prepareRelaxSpace();
    seedActinFilaments();
    seedMyosinMotors();
    /*ActinFilament* a = createNewFilament();
    a->direction[0] = 1;
    actins.push_back(a);
    ActinFilament* b = createNewFilament();
    b->direction[1] = 1;
    b->position[1] = 0.5*a->length + 0.25*a->diameter;
    actins.push_back(b);
    */
}

double Simulation::prepareForces(){

    for(MyosinMotorBinding* binding: bindings){

        binding->applyForces();

    }

    int N = actins.size();

    for(int i = 0; i<N; i++){
        Rod *rod = actins[i];
        for(int j = i+1; j<N; j++){
            Rod *other = actins[j];
            rod->collide(*other);
        }

        for(int j = 0; j<Constants::MYOSINS; j++){
            Rod *other = myosins[j];
            rod->collide(*other);
        }
    }

    for(int i = 0; i<Constants::MYOSINS; i++){
        Rod *rod = myosins[i];
        for(int j = i+1; j<Constants::MYOSINS; j++){
            Rod *other = myosins[j];
            rod->collide(*other);
        }
    }

    double sum = 0;
    for(int i = 0; i<N; i++){
        Rod *rod = actins[i];
        sum += rod->prepareForces();
    }

    for(int i = 0; i<Constants::MYOSINS; i++){
        Rod *rod = myosins[i];
        sum += rod->prepareForces();
    }

    return sum;
}

void Simulation::copyPositions(int index){
    int i = 0;
    std::vector<double> *positions = position_record[index].get();
    for(int a = 0; a<actins.size(); a++){
        Rod* rod = actins[a];
        for(int j = 0; j<3; j++){
            positions[0][i + j]= rod->position[j];
            positions[0][i + j + 3] = rod->direction[j];
        }
        i+=6;
    }

    for(MyosinMotor* rod : myosins){
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
    for(ActinFilament* rod : actins){
        for(int j = 0; j<3; j++){
            rod->position[j] = (*positions)[i + j];
            rod->direction[j] = (*positions)[i + j + 3];
        }
        rod->updateBounds();
        i+=6;
    }

    for(MyosinMotor* rod : myosins){
        for(int j = 0; j<3; j++){
            rod->position[j] = (*positions)[i + j];
            rod->direction[j] = (*positions)[i + j + 3];
        }
        rod->updateBounds();
        i+=6;
    }
}

void Simulation::copyForces(int index){
    int i = 0;
    std::vector<double> *forces = force_record[index].get();
    for(ActinFilament* rod : actins){
        for(int j = 0; j<3; j++){
            forces[0][i + j]= rod->force[j];
            forces[0][i + j + 3] = rod->torque[j];
        }
        i+=6;
    }

    for(MyosinMotor* rod : myosins){
        for(int j = 0; j<3; j++){
            forces[0][i + j] = rod->force[j];
            forces[0][i + j + 3] = rod->torque[j];
        }
        i+=6;
    }
}

void Simulation::clearForces(){

    for(ActinFilament* rod : actins){
        rod->clearForces();
    }

    for(MyosinMotor* rod : myosins){
        rod->clearForces();
    }
}

void Simulation::prepareForUpdate(int con_count, const std::vector<double> &coefficients) {
    int i = 0;
    for(ActinFilament* rod : actins){

        for(int k = 0; k<con_count; k++){

            double kn = coefficients[k];
            if(kn==0){continue;}

            std::vector<double> *forces = force_record[k].get();
            for(int j = 0; j<3; j++) {
                rod->force[j] += kn * forces[0][i + j];
                rod->torque[j] += kn * forces[0][i + j + 3];
            }
        }
        i+=6;
    }

    for(MyosinMotor* rod : myosins){
        for(int k = 0; k<con_count; k++){
            double kn = coefficients[k];
            if(kn==0){continue;}
            std::vector<double> *forces = force_record[k].get();
            for(int j = 0; j<3; j++) {
                rod->force[j] += kn * forces[0][i + j];
                rod->torque[j] += kn * forces[0][i + j + 3];
            }
        }

        i+=6;
    }
}

void Simulation::partialUpdate(double dt){
    for(ActinFilament* rod : actins){
        rod->update(dt);
        rod->clearForces();
        rod->updateBounds();
    }

    for(MyosinMotor* rod : myosins){
        rod->update(dt);
        rod->clearForces();
        rod->updateBounds();
    }
}

/**
* calculates root square of the difference between the current state (z_(k+1) and the more crude step, y_(k+1)
*/
double Simulation::calculateError(){
    int i = 0;
    double sum = 0;

    auto positions = position_record[1].get();
    double v;
    for(ActinFilament* rod : actins){
        for(int j = 0; j<3; j++){
            v = (*positions)[i + j] - rod->position[j];
            sum += v*v;
            v = (*positions)[i + j + 3] - rod->direction[j];
            sum+= v*v;
        }
        i+=6;
    }

    for(MyosinMotor* rod : myosins){
        for(int j = 0; j<3; j++){
            v = (*positions)[i + j] - rod->position[j];
            sum += v*v;
            v = (*positions)[i + j + 3] - rod->direction[j];
            sum+= v*v;
        }
        i+=6;
    }

    return sqrt(sum);
}
/**
* Using the Runge-Kutta-Fehlberg Method RKF45 an adaptive timestep technique
*
*/
void Simulation::relax(){
    bool relaxed = false;

    ActinFilament* a = actins[0];
    ActinFilament* b = actins[1];
    int stepped = 0;
    while(!relaxed) {
        double err, out_of_eq;
        do {
            //copy yk into position and direction record
            copyPositions(0);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //create k1/h
            out_of_eq = prepareForces();
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            copyForces(0);
            clearForces();
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //set the forces to be 1/(4h) k1
            prepareForUpdate(1, {0.25});
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //update the position to y0 + 1/4k1
            partialUpdate(working_dt);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //prepare k2
            prepareForces();
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            copyForces(1);
            clearForces();
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //move back to y0
            restorePositions(0);
            prepareForUpdate(2, {3.0 / 32.0, 9.0 / 32.0});
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //y = y0 + 3/32k1 + 9/32k2
            partialUpdate(working_dt);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //k3
            prepareForces();
            copyForces(2);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);
            clearForces();
            restorePositions(0);

            prepareForUpdate(3, {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0});
            //y = y0 + 1932/2197 k1 - 7200/2197 k2 + 7296/2197 k3
            partialUpdate(working_dt);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //k4
            prepareForces();
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);
            copyForces(3);
            clearForces();
            restorePositions(0);

            prepareForUpdate(4, {439.0 / 216.0, -8, 3680.0 / 513.0, -845.0 / 4104.0});
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //y = y0 + 439/216 k1 - 8 k2 + 3680/513 k3 - 845/4104 k4
            partialUpdate(working_dt);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //k5
            prepareForces();
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            copyForces(4);
            clearForces();
            restorePositions(0);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            prepareForUpdate(5, {-8.0 / 27.0, 2, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0});
            //y =y0 -8/27 k1 + 2 k2 - 3544/2565 k3 + 1859/4104 k4 - 11/40 k5
            partialUpdate(working_dt);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //k6 prepareForces();
            prepareForces();
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            copyForces(5);
            clearForces();
            restorePositions(0);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            prepareForUpdate(5, {25.0 / 216.0, 0, 1408.0 / 2565.0, 2197.0 / 4101.0, -0.2});
            partialUpdate(working_dt);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            //store y_(k+1) in the 2nd position.
            copyPositions(1);
            clearForces();
            restorePositions(0);


            prepareForUpdate(6, {16.0 / 135.0, 0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0});
            partialUpdate(working_dt);
            //printf("a: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", a->position[0], a->position[1], a->position[2], a->force[0], a->force[1], a->force[2]);
            //printf("b: %1.3e, %1.3e, %1.3e === %1.3e, %1.3e, %1.3e\n", b->position[0], b->position[1], b->position[2], b->force[0], b->force[1], b->force[2]);

            clearForces();

            err = calculateError();

            if(err>Constants::ERROR_THRESHOLD){
                restorePositions(0);
            }

            double v = (Constants::ERROR_THRESHOLD * working_dt*0.5/err);
            double new_dt = pow(v, 0.25)*working_dt;


            working_dt = new_dt>100*working_dt?100*working_dt:new_dt;
            working_dt = working_dt>Constants::DT?Constants::DT:working_dt;

        } while (err>Constants::ERROR_THRESHOLD);
        //printf("out of eq: %e, error: %e , dt %e\n", out_of_eq, err, working_dt);
        relaxed = out_of_eq<Constants::RELAXATION_LIMIT;
        stepped++;
        if(stepped>Constants::SUB_STEPS){
            printf("break\n");
            break;
        }
    }
    //printf("relaxed!\n");
}



void Simulation::step(){

    relax();
    for(ActinFilament* rod : actins){
        rod->updateBounds();
    }

    for(MyosinMotor* rod : myosins){
        rod->updateBounds();
    }

}

std::vector<ActinFilament*> &Simulation::getActins(){
    return actins;
}

std::vector<MyosinMotor*> &Simulation::getMyosins(){
    return myosins;
}