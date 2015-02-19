#include "Representations.h"

SpringRepresentation::SpringRepresentation(){
    //      floats per point x triangles x normals
    floats = 3*(2*rings*subdivisions + 2)*2*2;
}

void SpringRepresentation::setMaxSpringCount(int max){
    normal_offset = floats*max/2;
}

void SpringRepresentation::updateRepresentation(int index, float *positions, glm::dvec3 &a, glm::dvec3 &b) {
    int offset = index*floats/2;
    //printf("updating spring: %d\n", index);
    //printf("%f, %f, %f \t %f, %f, %f \n", a[0], a[1], a[2], b[0], b[1], b[2]);
    glm::dvec3 r = b - a;
    glm::dvec3 r_hat = glm::normalize(r);
    //printf("\tr^: %f,%f,%f\n", r_hat[0], r_hat[1], r_hat[2]);

    double min = 1;
    int axis = -1;

    for(int i = 0; i<3; i++){
        if(fabs(r_hat[i])<min){
            min = fabs(r_hat[i]);
            axis = i;
        }
    }

    glm::dvec3 j_hat(0,0,0);
    j_hat[axis] = 1;

    //printf("\tj^: %f,%f,%f\n", j_hat[0], j_hat[1], j_hat[2]);

    glm::dvec3 a1_hat = glm::normalize(glm::cross(r_hat, j_hat));
    glm::dvec3 a2_hat = glm::cross(a1_hat, r_hat);

    //printf("\ta1^: %f,%f,%f\n", a1_hat[0], a1_hat[1], a1_hat[2]);
    //printf("\ta2^: %f,%f,%f\n", a2_hat[0], a2_hat[1], a2_hat[2]);

    positions[offset + 0] = (float)a[0];
    positions[offset + 1] = (float)a[1];
    positions[offset + 2] = (float)a[2];
    positions[offset + 3] = (float)(a[0] + r_hat[0]*line_width);
    positions[offset + 4] = (float)(a[1] + r_hat[1]*line_width);
    positions[offset + 5] = (float)(a[2] + r_hat[2]*line_width);

    positions[normal_offset + offset + 0] = (float)a1_hat[0];
    positions[normal_offset + offset + 1] = (float)a1_hat[1];
    positions[normal_offset + offset + 2] = (float)a1_hat[2];
    positions[normal_offset + offset + 3] = (float)(a1_hat[0]);
    positions[normal_offset + offset + 4] = (float)(a1_hat[1]);
    positions[normal_offset + offset + 5] = (float)(a1_hat[2]);

    int n = rings*subdivisions;
    double ds = 1.0/n;
    double dtheta = 2*PI/subdivisions;
    for(int i = 0; i<2*n; i++){
        int stop = i/2;
        double theta = dtheta*stop;
        double s = ds*stop;

        positions[offset + 6*(i+1)] = (float)(a[0] + r[0]*s + a1_hat[0]*radius*cos(theta) + a2_hat[0]*radius*sin(theta));
        positions[offset + 6*(i+1) + 1] = (float)(a[1] + r[1]*s + a1_hat[1]*radius*cos(theta) + a2_hat[1]*radius*sin(theta));
        positions[offset + 6*(i+1) + 2] = (float)(a[2] + r[2]*s + a1_hat[2]*radius*cos(theta) + a2_hat[2]*radius*sin(theta));
        positions[offset + 6*(i+1)+3] = positions[offset + 6*(i+1)] + (float)(r_hat[0]*line_width);
        positions[offset + 6*(i+1)+4] = positions[offset + 6*(i+1) + 1] + (float)(r_hat[1]*line_width);
        positions[offset + 6*(i+1)+5] = positions[offset + 6*(i+1) + 2] + (float)(r_hat[2]*line_width);

        positions[normal_offset + offset + 6*(i+1)] = (float)(a1_hat[0]*cos(theta) + a2_hat[0]*sin(theta));
        positions[normal_offset + offset + 6*(i+1) + 1] = (float)(a1_hat[1]*cos(theta) + a2_hat[1]*sin(theta));
        positions[normal_offset + offset + 6*(i+1) + 2] = (float)(a1_hat[2]*cos(theta) + a2_hat[2]*sin(theta));
        positions[normal_offset + offset + 6*(i+1)+3] = positions[normal_offset + offset + 6*(i+1)];
        positions[normal_offset + offset + 6*(i+1)+4] = positions[normal_offset + offset + 6*(i+1) + 1];
        positions[normal_offset + offset + 6*(i+1)+5] = positions[normal_offset + offset + 6*(i+1) + 2];

    }

    positions[offset + 6*(2*n + 1) + 0] = (float)b[0];
    positions[offset + 6*(2*n + 1) + 1] = (float)b[1];
    positions[offset + 6*(2*n + 1) + 2] = (float)b[2];


    positions[offset + 6*(2*n + 1) + 3] = (float)b[0]+ (float)(r_hat[0]*line_width);
    positions[offset + 6*(2*n + 1) + 4] = (float)b[1] + (float)(r_hat[1]*line_width);
    positions[offset + 6*(2*n + 1) + 5] = (float)b[2] + (float)(r_hat[2]*line_width);

    int first = normal_offset + offset + 6*(2*n + 1);
    double theta = 0;
    positions[first] = (float)(a1_hat[0]*cos(theta) + a2_hat[0]*sin(theta));
    positions[first + 1] = (float)(a1_hat[1]*cos(theta) + a2_hat[1]*sin(theta));
    positions[first + 2] = (float)(a1_hat[2]*cos(theta) + a2_hat[2]*sin(theta));
    positions[first+3] = positions[first];
    positions[first+4] = positions[first + 1];
    positions[first+5] = positions[first + 2];
    /*
    printf("floats: %d \t\t %d\n", floats, first+6);
    float* t=positions;
    for(int i = 0; i<floats/6; i++){

        printf("%d::\t %f, %f, %f\n",i, t[offset + 3*i+0], t[offset + 3*i+1], t[offset + 3*i+2],
                t[normal_offset + offset + 3*i+0], t[normal_offset + offset + 3*i+1], t[normal_offset + offset + 3*i+2]);

    }
    */
}



int SpringRepresentation::getFloatCount(){
    return floats;
}
