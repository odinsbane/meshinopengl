#include "Representations.h"

int RectangularPrism::getFloatCount(){
    //N boxes with 6 sides, 2 triangles per side, 3 nodes per triangle, 3 position and 3 normal per node.
    // N          *6       *2                    *3                   * (3              +3) =
    return N*6*2*3*(3 + 3);
}

int RectangularPrism::getPositionOffset(){
    return N*SIDES*TRIANGLES*NODES*POSITIONS;
}

void RectangularPrism::updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c){
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

    glm::dvec3 norm = glm::cross(b - a, c - a);
    norm = glm::normalize(norm);
    target[0 + position_offset] = static_cast<float>(norm[0]);
    target[1 + position_offset] = static_cast<float>(norm[1]);
    target[2 + position_offset] = static_cast<float>(norm[2]);
    target[3 + position_offset] = static_cast<float>(norm[0]);
    target[4 + position_offset] = static_cast<float>(norm[1]);
    target[5 + position_offset] = static_cast<float>(norm[2]);
    target[6 + position_offset] = static_cast<float>(norm[0]);
    target[7 + position_offset] = static_cast<float>(norm[1]);
    target[8 + position_offset] = static_cast<float>(norm[2]);

}

void RectangularPrism::updateFace(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &d){
    int HALF = NODES*POSITIONS;
    updateTriangle(target, a, b, c);
    updateTriangle(&target[HALF], a, c, d);

}

int RectangularPrism::getElementNodeCount(){
    return NODES*TRIANGLES*SIDES;
}

void RectangularPrism::updateRod(int start, Rod &rod){
    int FACE = TRIANGLES*NODES*POSITIONS;

    //index times the number of floats for a rod.
    //int start = index*SIDES*FACE;
    //step one create the 3 axis for changes
    glm::dvec3 main_axis(
            rod.direction[0]*0.5*(rod.length + rod.diameter),
            rod.direction[1]*0.5*(rod.length + rod.diameter),
            rod.direction[2]*0.5*(rod.length + rod.diameter)
    );

    glm::dvec3 first_axis = glm::cross(main_axis, zaxis);
    if(glm::length(first_axis)==0){
        first_axis = glm::cross(main_axis, xaxis);
    }
    first_axis = glm::normalize(first_axis);
    first_axis[0] = first_axis[0]*0.5*rod.diameter;
    first_axis[1] = first_axis[1]*0.5*rod.diameter;
    first_axis[2] = first_axis[2]*0.5*rod.diameter;

    glm::dvec3 second_axis = glm::cross(main_axis, first_axis);
    second_axis = glm::normalize(second_axis);
    second_axis[0] = second_axis[0]*0.5*rod.diameter;
    second_axis[1] = second_axis[1]*0.5*rod.diameter;
    second_axis[2] = second_axis[2]*0.5*rod.diameter;
    {
        //front
        glm::dvec3 a = rod.position + main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position + main_axis - first_axis + second_axis;
        glm::dvec3 c = rod.position + main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position + main_axis + first_axis - second_axis;

        updateFace(&positions[start], a, b, c, d);
    }

    {
        //back
        glm::dvec3 a = rod.position - main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position - main_axis + first_axis - second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position - main_axis - first_axis + second_axis;

        updateFace(&positions[start + FACE], a, b, c, d);
    }

    {
        //right
        glm::dvec3 a = rod.position + main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position + main_axis + first_axis - second_axis;
        glm::dvec3 c = rod.position - main_axis + first_axis - second_axis;
        glm::dvec3 d = rod.position - main_axis + first_axis + second_axis;

        updateFace(&positions[start + 2*FACE], a, b, c, d);
    }

    {
        //left
        glm::dvec3 a = rod.position + main_axis - first_axis + second_axis;
        glm::dvec3 b = rod.position - main_axis - first_axis + second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position + main_axis - first_axis - second_axis;

        updateFace(&positions[start + 3*FACE], a, b, c, d);
    }

    {
        //top
        glm::dvec3 a = rod.position + main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position - main_axis + first_axis + second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis + second_axis;
        glm::dvec3 d = rod.position + main_axis - first_axis + second_axis;

        updateFace(&positions[start + 4*FACE], a, b, c, d);
    }

    {
        //bottom
        glm::dvec3 a = rod.position + main_axis + first_axis - second_axis;
        glm::dvec3 b = rod.position + main_axis - first_axis - second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position - main_axis + first_axis - second_axis;

        updateFace(&positions[start + 5*FACE], a, b, c, d);
    }

    /*printf("box \n\n");
    for(int i = 0; i<36; i++){
        float* p = &positions[3*i];
        float* n = &positions[3*i + position_offset];
        printf("%1.1f, %1.1f, %1.1f, ... %1.1f, %1.1f, %1.1f \n", p[0], p[1], p[2], n[0], n[1], n[2]);
    }*/
}
