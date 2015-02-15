#include "Representations.h"

MeshMyosin::MeshMyosin(int rods) {
    // a panel is a 'square' shape composed of two triangles.
    //        * stalk c_div x l_divs panels.
    //        * top&bottom c_div panels + c_div triangles
    element_node_count = 6*c_divs*l_divs + 2*(c_divs*6 + c_divs*3);
    floats = rods*element_node_count*6; //3 positions per node, 3 normals per node
    position_offset = floats/2;
}

int MeshMyosin::getFloatCount(){
    return floats;
}

const glm::dmat3x3 identity(1,0,0,0,1,0,0,0,1);

void MeshMyosin::updateRod(int start, Rod &rod){
    int p_triangle = 3*3;
    int p_face = 2*p_triangle;
    int p_strip = c_divs*p_face;
    //create our transform.
    glm::dmat3x3 rotator;
    double dot_p = glm::dot(rod.direction, zaxis);
    if(dot_p == 1){
        rotator = identity;
    } else if(dot_p==-1){
        rotator = identity;
        rotator[2][2] = -1;
    } else {
        glm::dvec3 rot_axis = glm::cross(zaxis, rod.direction);
        rot_axis = glm::normalize(rot_axis);
        double angle = acos(dot_p);
        glm::dmat4x4 more = glm::rotate(angle, rot_axis);
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                rotator[i][j] = more[i][j];
            }
        }
    }
    /*
    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){

            printf("\t%2.1f", rotator[i][j] );

        }
        printf("\n");
    }
    */

    //stalk
    glm::dvec3 a,b,c,d,na,nb,nc,nd;
    double dtheta = 2*PI/c_divs;
    double dl = rod.length/l_divs;
    double theta;
    double ntheta;
    double radius = rod.diameter*0.5;
    double length = rod.length;

    //the stalk
    for(int j = 0; j<l_divs; j++) {
        double s = j*dl;
        for (int i = 0; i < c_divs; i++) {
            theta = i*dtheta;
            ntheta = ((i+1)%c_divs) * (dtheta);
            a = rotator*getStalkPosition(s, theta, radius, length) + rod.position;
            b = rotator*getStalkPosition(s+dl, theta, radius, length) + rod.position;
            c = rotator*getStalkPosition(s+dl, ntheta, radius, length) + rod.position;
            d = rotator*getStalkPosition(s, ntheta, radius, length) + rod.position;

            na = rotator*getStalkNormal(s, theta, radius, length);
            nb = rotator*getStalkNormal(s+dl, theta, radius, length);
            nc = rotator*getStalkNormal(s+dl, ntheta, radius, length);
            nd = rotator*getStalkNormal(s, ntheta, radius, length);
            /*
            printf("%f\t%f\t%f\t", a[0], a[1], a[2]);
            printf("%f\t%f\t%f\t", b[0], b[1], b[2]);
            printf("%f\t%f\t%f\t", c[0], c[1], c[2]);
            printf("%f\t%f\t%f\n", d[0], d[1], d[2]);
            */
            updateTriangle(&positions[start + i*p_face + j*p_strip], b, a, d, nb, na, nd);
            updateTriangle(&positions[start + i*p_face + j*p_strip + p_triangle], c, b, d, nc, nb, nd);
        }
    }

    start = start + l_divs*p_strip;

    for(int i = 0; i<c_divs; i++){
        double s = -0.5*radius;
        theta = i*dtheta;
        ntheta = ((i+1)%c_divs) * (dtheta);
        a = rotator*getTipPosition(s, theta, radius, length) + rod.position;
        b = rotator*getTipPosition(0, theta, radius, length) + rod.position;
        c = rotator*getTipPosition(0, ntheta, radius, length) + rod.position;
        d = rotator*getTipPosition(s, ntheta, radius, length) + rod.position;

        na = rotator*getTipNormal(s, theta, radius, length);
        nb = rotator*getTipNormal(0, theta, radius, length);
        nc = rotator*getTipNormal(0, ntheta, radius, length);
        nd = rotator*getTipNormal(s, ntheta, radius, length);

        updateTriangle(&positions[start + i*p_face], b, a, d, nb, na, nd);
        updateTriangle(&positions[start + i*p_face + p_triangle], c, b, d, nc, nb, nd);
    }

    start = start + p_strip;

    for(int i = 0; i<c_divs; i++){
        double s = -radius;
        theta = i*dtheta;
        ntheta = ((i+1)%c_divs) * (dtheta);
        b = rotator*getTipPosition(0.5*s, theta, radius, length) + rod.position;
        c = rotator*getTipPosition(0.5*s, ntheta, radius, length) + rod.position;
        d = rotator*getTipPosition(s, ntheta, radius, length) + rod.position;

        nb = rotator*getTipNormal(0.5*s, theta, radius, length);
        nc = rotator*getTipNormal(0.5*s, ntheta, radius, length);
        nd = rotator*getTipNormal(s, ntheta, radius, length);


        updateTriangle(&positions[start + i*p_triangle], c, b, d, nc, nb, nd);
    }

    start = start + p_strip/2;

    for(int i = 0; i<c_divs; i++){
        double s = length;
        theta = i*dtheta;
        ntheta = ((i+1)%c_divs) * (dtheta);
        a = rotator*getTipPosition(s, theta, radius, length) + rod.position;
        b = rotator*getTipPosition(s+radius*0.5, theta, radius, length) + rod.position;
        c = rotator*getTipPosition(s+radius*0.5, ntheta, radius, length) + rod.position;
        d = rotator*getTipPosition(s, ntheta, radius, length) + rod.position;

        na = rotator*getTipNormal(s, theta, radius, length);
        nb = rotator*getTipNormal(s+0.5*radius, theta, radius, length);
        nc = rotator*getTipNormal(s+0.5*radius, ntheta, radius, length);
        nd = rotator*getTipNormal(s, ntheta, radius, length);


        updateTriangle(&positions[start + i*p_face], b, a, d, nb, na, nd);
        updateTriangle(&positions[start + i*p_face + p_triangle], c, b, d, nc, nb, nd);
    }

    start = start + p_strip;

    for(int i = 0; i<c_divs; i++){
        double s = length + 0.5*radius;
        theta = i*dtheta;
        ntheta = ((i+1)%c_divs) * (dtheta);
        a = rotator*getTipPosition(s, theta, radius, length) + rod.position;
        b = rotator*getTipPosition(s+radius*0.5, theta, radius, length) + rod.position;
        d = rotator*getTipPosition(s, ntheta, radius, length) + rod.position;

        na = rotator*getTipNormal(s, theta, radius, length);
        nb = rotator*getTipNormal(s+0.5*radius, theta, radius, length);
        nd = rotator*getTipNormal(s, ntheta, radius, length);


        updateTriangle(&positions[start + i*p_triangle], b, a, d, nb, na, nd);

    }
}

int MeshMyosin::getPositionOffset(){
    return position_offset;
}
int MeshMyosin::getElementNodeCount(){
    return element_node_count;
}

glm::dvec3 MeshMyosin::getStalkPosition(double s, double theta, double radius, double length){
    double alpha = pitch*s;
    double a = (1+eccentricity)*radius*cos(theta - alpha);
    double b = (radius)*sin(theta-alpha);
    double r = sqrt(a*a + b*b);
    double y = sin(theta)*r;
    double x = cos(theta)*r;
    double z = s - length/2.0;
    return glm::dvec3(x,y,z);

}

glm::dvec3 MeshMyosin::getTipPosition(double s, double theta, double radius, double length){
    double alpha = s>0?pitch*length:0;

    double a = (1+eccentricity)*radius*cos(theta - alpha);
    double b = (radius)*sin(theta-alpha);
    double r = sqrt(a*a + b*b);
    if(s<=0){
        if(-s>0.9*radius){
            r = 0;
        } else{

            r = sqrt(r*r - s*s);
        }
    } else{
        double h = s - length;
        if(h>0.9*radius){
            r = 0;
        } else{

            r = sqrt(r*r - h*h);
        }

    }

    double y = sin(theta)*r;
    double x = cos(theta)*r;
    double z = s - length/2.0;
    return glm::dvec3(x,y,z);
}
const double sqrt2 = sqrt(2)/2.0;
glm::dvec3 MeshMyosin::getTipNormal(double s, double theta, double radius, double length){
    double alpha = s>0?pitch*length:0;

    double a = (1+eccentricity)*radius*cos(theta - alpha);
    double b = (radius)*sin(theta-alpha);
    double r = sqrt(a*a + b*b);

    double ap = -(eccentricity + 1)*radius*sin(theta - alpha);
    double bp = radius*cos(theta-alpha);

    double rprime = (a*ap + b*bp)/r;

    double xp = -r*sin(theta) + cos(theta)*rprime;
    double yp = r*cos(theta) + sin(theta)*rprime;
    double n = sqrt(xp*xp + yp*yp);

    double z = 0;

    if(s<=0){
        if(-s>0.9*radius){
            xp = 0;
            yp = 0;
            z = -1;
        } else{
            z = s/r;

            xp = xp*(r+s)/r;
            yp = yp*(r+s)/r;
        }
    } else{
        if(s-length>0.9*radius){
            xp = 0;
            yp = 0;
            z = 1;
        } else{
            double h = s -length;
            z = h/r;

            xp = xp*(r-h)/r;
            yp = yp*(r-h)/r;

        }

    }

    double x = yp/n;
    double y = -xp/n;
    return glm::dvec3(x,y,z);
}

glm::dvec3 MeshMyosin::getStalkNormal(double s, double theta, double radius, double length){
    double alpha = pitch*s;
    //double a = (eccentricity + 1)*radius)*cos(theta - alpha);
    //double b = (radius)*sin(theta-alpha);
    double a = (1+eccentricity)*radius*cos(theta - alpha);
    double b = (radius)*sin(theta-alpha);
    double r = sqrt(a*a + b*b);

    double ap = -(eccentricity + 1)*radius*sin(theta - alpha);
    double bp = radius*cos(theta-alpha);

    double rprime = (a*ap + b*bp)/r;
    double xp = -r*sin(theta) + cos(theta)*rprime;
    double yp = r*cos(theta) + sin(theta)*rprime;
    double n = sqrt(xp*xp + yp*yp);
    double x = yp/n;
    double y = -xp/n;
    double z = 0;
    return glm::dvec3(x,y,z);
}

void MeshMyosin::updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &na, glm::dvec3 &nb, glm::dvec3 &nc){
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

void MeshMyosin::setPositionOffset(int offset){
    position_offset = offset;
}
