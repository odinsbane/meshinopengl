#include "interactions.h"

glm::dvec3 getReflectedPoint(glm::dvec3 &src, glm::dvec3 &target) {
    glm::dvec3 out;
    double hw = Constants::WIDTH*0.5;
    if( target[0] - src[0] > hw){
        out[0] = target[0] - Constants::WIDTH;
    } else if(target[0] - src[0]<-hw){
        out[0] = target[0] + Constants::WIDTH;
    } else{
        out[0] = target[0];
    }


    if( target[1] - src[1] > hw){
        out[1] = target[1] - Constants::WIDTH;
    } else if(target[1] - src[1]<-hw){
        out[1] = target[1] + Constants::WIDTH;
    } else{
        out[1] = target[1];
    }

    out[2] = target[2];
    return out;
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

void MyosinMotorBinding::update(double dt){
    checkTime(MyosinMotor::FRONT, dt);
    checkTime(MyosinMotor::BACK, dt);

    slideHead(MyosinMotor::FRONT, dt);
    slideHead(MyosinMotor::BACK, dt);


    if(Constants::MYOSIN_DIFFUSION_FORCE>0){

        motor->clearForces();

        double f = 2*Constants::MYOSIN_DIFFUSION_FORCE/(sqrt(3));

        motor->applyForce( new glm::dvec4(
                f*(number_generator->nextDouble() - 0.5),
                f*(number_generator->nextDouble() - 0.5),
                f*(number_generator->nextDouble() - 0.5),
                motor->length*(number_generator->nextDouble() - 0.5)  )
        );
        motor->prepareForces();
        motor->update(dt);
    }


    if(
            fabs(motor->getPoint(motor->length*0.5)[2])>Constants::THICKNESS*0.5 &&
                    fabs(motor->getPoint(-motor->length*0.5)[2])>Constants::THICKNESS*0.5
            ){

        motor->unbind(MyosinMotor::FRONT);
        motor->unbind(MyosinMotor::BACK);
        //model->placeBoundMyosinMotor(motor, this);
    }
}

void MyosinMotorBinding::checkTime(int head, double dt){
    if(motor->isBound(head)){
        current_time[head] += dt;
        if(current_time[head]>unbind_time[head]){
            motor->unbind(head);
        }
    }
}

void MyosinMotorBinding::slideHead(int head, double dt) {
    if(motor->isFree(head)) return;
    binding_position[head] += sliding[head]*dt;

    if(fabs(binding_position[head])>motor->getBound(head)->length*0.5) {
        //unbind
        motor->unbind(head);
    }
    //printf("slide: %16.16e \n", sliding[head]*dt);
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
    finito = true;
    filaments[0]->unbind(filaments[1]);
    filaments[1]->unbind(filaments[0]);
}

bool CrosslinkedFilaments::finished(){
    return finito;
}

glm::dvec3 MyosinMotorBinding::getBindingPosition(int head) {
    if(motor->isBound(head)){
        return motor->getBound(head)->getPoint(binding_position[head]);
    } else{
        double s = 0.5*motor->length - Constants::MYOSIN_BIND_LENGTH;
        if(head==MyosinMotor::BACK) s = -s;
        return motor->getPoint(s);
    }


}

glm::dvec3 MyosinMotorBinding::getHeadPosition(int head) {
    double s = 0.5*motor->length;
    if(head==MyosinMotor::BACK) s = -s;
    return motor->getPoint(s);
}
