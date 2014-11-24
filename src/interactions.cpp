#include "interactions.h"

glm::dvec3 getReflectedPoint(glm::dvec3 &src, glm::dvec3 &target){
    return target;
}
void MyosinMotorBinding::bind(ActinFilament* f, int head, double position){
    motor->bind(f, head);
    binding_position[head] = position;
    current_time[head] = 0;
    unbind_time[head] = -motor->tau_B*log(number_generator->nextDouble());
}

void MyosinMotorBinding::applyForces(){
    if(motor->isBound(MyosinMotor::FRONT)){
        headForce(MyosinMotor::FRONT);
    }

    if(motor->isBound(MyosinMotor::BACK)){
        headForce(MyosinMotor::BACK);
    }

}

void MyosinMotorBinding::headForce(int head){
    ActinFilament* filament = motor->getBound(head);
    double ml = head==MyosinMotor::FRONT?0.5*motor->length:-0.5*motor->length;
    glm::dvec3 front_head = motor->getPoint(ml);
    glm::dvec3 attached_position = filament->getPoint( binding_position[head]);

    glm::dvec3 reflected = getReflectedPoint(front_head, attached_position);

    //line from myosin head to actin filament attachment location.
    glm::dvec3 r = reflected - front_head;

    double mag = glm::length(r);
    double separation = mag - motor->bind_length;
    //crude approximation needs to be fixed. TODO
    double f = mag==0?0:motor->K_m/mag*separation;

    glm::dvec4* head_force = new glm::dvec4( f*r[0], f*r[1], f*r[2], ml );
    motor->applyForce(head_force);

    double dot = (*head_force)[0]*filament->direction[0] + (*head_force)[1]*filament->direction[1] + (*head_force)[2]*filament->direction[2];
    sliding[head] = (motor->F0 - dot)/motor->alpha_s;

    glm::dvec4* filament_force = new glm::dvec4(-(*head_force)[0], -(*head_force)[1], -(*head_force)[2], binding_position[head]);
    filament->applyForce(filament_force);
}

void CrosslinkedFilaments::applyForces() {

    glm::dvec3 a_pos = filaments[0]->getPoint(locations[0]);
    glm::dvec3 b_pos = filaments[1]->getPoint(locations[1]);
    b_pos = getReflectedPoint(a_pos, b_pos);
    glm::dvec3 r = b_pos - a_pos;
    double mag = glm::length(r);
    double f = 0;

    if(mag>0){
        f = (mag - length)*K_x/mag;
    }
    glm::dvec4* force_a = new glm::dvec4(r[0]*f, r[1]*f, r[2]*f, locations[0]);
    filaments[0]->applyForce(force_a);
    glm::dvec4* force_b = new glm::dvec4(-r[0]*f, -r[1]*f, -r[2]*f, locations[1]);
    filaments[1]->applyForce(force_b);


}

void CrosslinkedFilaments::update(double dt){
    current_time += dt;
    if(current_time>unbind_time){
        unbind();
    }
}

void CrosslinkedFilaments::unbind() {

}

