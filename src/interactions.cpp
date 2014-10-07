#include "interactions.h"

void MyosinMotorBinding::bind(ActinFilament* f, int head, double position){
    motor->bind(f, head);
    binding_position[head] = position;
    current_time[head] = 0;
    unbind_time[head] = -motor->tau_B*log(number_generator->nextDouble());
}

void MyosinMotorBinding::applyForces(){


}