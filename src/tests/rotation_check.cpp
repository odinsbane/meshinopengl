#include <iostream>
#define GLM_FORCE_RADIANS
#include <glm/geometric.hpp> //glm::dot
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/gtx/transform.hpp>

using namespace std;
int main(int arg_c, char** args){
    glm::dvec3 rot(1,0,0);
    glm::dmat4x4 mat = glm::rotate(0.0, rot);
    glm::dmat3x3 om;
    for(int i = 0; i<4; i++) {
        printf(". %f %f %f %f\n", mat[0][i], mat[1][i], mat[2][i], mat[3][i]);
    }

    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            om[i][j] = mat[i][j];
        }
    }


    glm::dvec3 before(0,1,0);
    glm::dvec3 after = om*before;
    printf("result %f %f %f \n", after[0], after[1], after[2]);
    return 0;
}