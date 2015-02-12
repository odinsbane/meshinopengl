#include "Representations.h"

MeshHelix::MeshHelix(int rods) {
    // a panel is a 'square' shape composed of two triangles.
    //        * stalk c_div x l_divs panels.
    //        * top&bottom c_div panels + c_div triangles
    element_node_count = 6*c_divs*l_divs + 2*(c_divs*6 + c_divs*3);
    floats = rods*element_node_count*6; //3 positions per node, 3 normals per node
    position_offset = floats/2;
}

int MeshHelix::getFloatCount(){
    return floats;
}

const glm::dmat3x3 identity(1,0,0,0,1,0,0,0,1);

void MeshHelix::updateRod(int index, Rod &rod){
    //create our transform.
    glm::dmat3x3 rotator;
    double dot_p = glm::dot(rod.direction, zaxis);

    if(dot_p == 1){
        rotator = identity;
    } else if(dot_p==-1){
        rotator = identity;
        rotator[2][2] = -1;
    } else {
        glm::dvec3 rot_axis = glm::cross(rod.direction, zaxis);
        rot_axis = glm::normalize(rot_axis);
        double angle = acos(dot_p);
        glm::dmat4x4 more = glm::rotate(angle, rot_axis);
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                rotator[i][j] = more[i][j];
            }
        }
    }


    //stalk
    glm::dvec3 a,b,c,d,na,nb,nc,nd;
    double dtheta = 2*PI/c_divs;
    double dl = rod->length/l_divs;
    for(int j = 0; j<l_divs; j++) {
        double s = j*dl;
        for (int i = 0; i < c_divs; i++) {

        }
    }

}

int MeshHelix::getPositionOffset(){
    return position_offset;
}
int MeshHelix::getElementNodeCount(){
    return element_node_count;
}

void MeshHelix::getStatus(glm::dvec3& pos, glm::dvec3& norm, double s, double theta, double radius){

}