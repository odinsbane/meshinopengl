#include "Representations.h"

MeshCylinder::MeshCylinder(int rods, int divisions){
    //        stalk points.         conical caps
    floats = rods*divisions*2*3*6 + 2*rods*divisions*3*(3+3);
    position_offset = floats/2;
    element_node_count = divisions*2*3 + 2*divisions*3;
    this->divisions = divisions;
}

int MeshCylinder::getFloatCount(){
    return floats;
}

int MeshCylinder::getPositionOffset(){
    return position_offset;
}

int MeshCylinder::getElementNodeCount() {
    return element_node_count;
}

void MeshCylinder::updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &na, glm::dvec3 &nb, glm::dvec3 &nc){
    target[0] = static_cast<float>(a[0]);
    target[1] = static_cast<float>(a[1]);
    target[2] = static_cast<float>(a[2]);
    target[3] = static_cast<float>(b[0]);
    target[4] = static_cast<float>(b[1]);
    target[5] = static_cast<float>(b[2]);
    target[6] = static_cast<float>(c[0]);
    target[7] = static_cast<float>(c[1]);
    target[8] = static_cast<float>(c[2]);

    int position_offset = getPositionOffset();

    target[0 + position_offset] = static_cast<float>(na[0]);
    target[1 + position_offset] = static_cast<float>(na[1]);
    target[2 + position_offset] = static_cast<float>(na[2]);
    target[3 + position_offset] = static_cast<float>(nb[0]);
    target[4 + position_offset] = static_cast<float>(nb[1]);
    target[5 + position_offset] = static_cast<float>(nb[2]);
    target[6 + position_offset] = static_cast<float>(nc[0]);
    target[7 + position_offset] = static_cast<float>(nc[1]);
    target[8 + position_offset] = static_cast<float>(nc[2]);

}
glm::dvec3 sum(double a, glm::dvec3 &va, double b, glm::dvec3 &vb){

    return glm::dvec3(
            a*va[0] + b*vb[0],
            a*va[1] + b*vb[1],
            a*va[2] + b*vb[2]
    );


}
void MeshCylinder::updateRod(int index, Rod &rod) {

    //index times the number of floats for a rod.
    int start = index*element_node_count*3;
    //step one create the 3 axis for changes
    glm::dvec3 main_axis(
            rod.direction[0]*0.5*(rod.length),
            rod.direction[1]*0.5*(rod.length),
            rod.direction[2]*0.5*(rod.length)
    );

    glm::dvec3 first_axis = glm::cross(main_axis, zaxis);
    if(glm::length(first_axis)==0){
        first_axis = glm::cross(main_axis, xaxis);
    }
    first_axis = glm::normalize(first_axis);

    glm::dvec3 second_axis = glm::cross(main_axis, first_axis);
    second_axis = glm::normalize(second_axis);

    double dtheta = 2*3.14159/divisions;

    double r = 0.5* rod.diameter;
    int floats_per_face = 3*6 + 6*3;

    glm::dvec3 tip = sum(0.5*rod.diameter + 0.5*rod.length, rod.direction, 1, rod.position);
    glm::dvec3 tail = sum(-0.5*(rod.diameter + rod.length), rod.direction, 1, rod.position);
    glm::dvec3 backwards(-rod.direction[0], -rod.direction[1], -rod.direction[2]);
    for(int i = 0;i<divisions; i++){
        double a1 = cos(dtheta*i);
        double b1 = sin(dtheta*i);
        double a2 = cos(dtheta*((i+1)%divisions));
        double b2 = sin(dtheta*((i+1)%divisions));

        glm::dvec3 n1 = sum(a1, first_axis, b1, second_axis);
        glm::dvec3 n2 = sum(a2, first_axis, b2, second_axis);

        glm::dvec3 a = sum(-1, main_axis, r, n1);
        glm::dvec3 b = sum(-1, main_axis, r, n2);
        glm::dvec3 c = sum(1, main_axis, r, n2);
        glm::dvec3 d = sum(1, main_axis, r, n1);

        a = sum(1, rod.position, 1, a);
        b = sum(1, rod.position, 1, b);
        c = sum(1, rod.position, 1, c);
        d = sum(1, rod.position, 1, d);

        updateTriangle(&positions[start + i*floats_per_face], a, b, d, n1, n2, n1);
        updateTriangle(&positions[start + i*floats_per_face + 9], b, c, d, n2, n2, n1);

        updateTriangle(&positions[start + i*floats_per_face + 18], d, c, tip, n1, n2, rod.direction);
        updateTriangle(&positions[start + i*floats_per_face + 27], b, a, tail, n2, n1, backwards);
    }

}
