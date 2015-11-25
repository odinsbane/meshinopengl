#ifndef __ParallelBalls__Representations__
#define __ParallelBalls__Representations__
#include "rod.h"
#include <glm/geometric.hpp> //glm::dot
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

const glm::dvec3 zaxis(0,0,1);
const glm::dvec3 xaxis(1,0,0);

class CylinderRepresentation{

    public:
    CylinderRepresentation(){}
    float* positions;
    virtual int getFloatCount(){ return 0;}
    virtual int getPositionOffset(){return 0;}
    virtual void setPositionOffset(int offset){}
    virtual void updateRod(int start, Rod &rod){}
    virtual int getElementNodeCount(){return 0;}
    void setPositions(float * p){positions = p;}
};

class RectangularPrism : public CylinderRepresentation{
    //divisions.
    const int SIDES=6;
    const int TRIANGLES=2;
    const int NODES=3;
    const int POSITIONS=3;
    const int NORMALS=3;
    int N;
    void updateFace(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &d);
    void updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c);
public:
    RectangularPrism(int rods){
        N=rods;
    }
    int getFloatCount();
    void updateRod(int start, Rod &rod);
    int getPositionOffset();
    int getElementNodeCount();
};

class MeshCylinder : public CylinderRepresentation{
    int divisions = 16;
    int floats, position_offset, element_node_count;
    void updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &na, glm::dvec3 &nb, glm::dvec3 &nc);
public:
    MeshCylinder(int rods);
    int getFloatCount();
    void updateRod(int start, Rod &rod);
    int getPositionOffset();
    int getElementNodeCount();
    void setPositionOffset(int offset);
};

class MeshHelix : public CylinderRepresentation{
    int c_divs = 5;
    int l_divs = 10;
    double pitch = 42.4539; //2 pi every 37nm
    //double pitch = PI/2.0;
    double eccentricity = 0.5;
    int floats, position_offset, element_node_count;
    void updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &na, glm::dvec3 &nb, glm::dvec3 &nc);
    public:
        MeshHelix(int rods);
        int getFloatCount();
        void updateRod(int start, Rod &rod);
        int getPositionOffset();
        void setPositionOffset(int offset);
        int getElementNodeCount();
        glm::dvec3 getStalkPosition(double s, double theta, double radius, double length);
        glm::dvec3 getStalkNormal(double s, double theta, double radius, double length);
        glm::dvec3 getTipPosition(double s, double theta, double radius, double length);
        glm::dvec3 getTipNormal(double s, double theta, double radius, double length);
};

class MeshMyosin : public CylinderRepresentation{
    int c_divs = 50;
    int l_divs = 20;
    double pitch = 0; //2 pi every 37nm
    //double pitch = PI/2.0;
    double eccentricity = 0.1;
    int floats, position_offset, element_node_count;
    void updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &na, glm::dvec3 &nb, glm::dvec3 &nc);
public:
    MeshMyosin(int rods);
    int getFloatCount();
    void updateRod(int start, Rod &rod);
    int getPositionOffset();
    void setPositionOffset(int offset);
    int getElementNodeCount();
    glm::dvec3 getStalkPosition(double s, double theta, double radius, double length);
    glm::dvec3 getStalkNormal(double s, double theta, double radius, double length);
    glm::dvec3 getTipPosition(double s, double theta, double radius, double length);
    glm::dvec3 getTipNormal(double s, double theta, double radius, double length);
};

class SpringRepresentation{
    int floats;
    int normal_offset;
    int rings = 4;
    int subdivisions = 28;
    double radius = 0.035;
    double line_width = 0.03;

public:
    SpringRepresentation();
    void updateRepresentation(int index, float* positions, glm::dvec3 &a, glm::dvec3 &b);
    int getFloatCount();
    void setMaxSpringCount(int number);
};

#endif
