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
    virtual void updateRod(int index, Rod &rod){}
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
    void updateRod(int index, Rod &rod);
    int getPositionOffset();
    int getElementNodeCount();
};

class MeshCylinder : public CylinderRepresentation{
    int divisions, floats, position_offset, element_node_count;
    void updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &na, glm::dvec3 &nb, glm::dvec3 &nc);
public:
    MeshCylinder(int rods, int divisions);
    int getFloatCount();
    void updateRod(int index, Rod &rod);
    int getPositionOffset();
    int getElementNodeCount();

};

class MeshHelix : public CylinderRepresentation{
    int c_divs = 20;
    int l_divs = 20;
    int floats, position_offset, element_node_count;
    public:
        MeshHelix(int rods);
        int getFloatCount();
        void updateRod(int index, Rod &rod);
        int getPositionOffset();
        int getElementNodeCount();
        void getStatus(glm::dvec3& pos, glm::dvec3& norm, double s, double theta, double radius);
};

class SpringRepresentation{
    int floats;
    int rings = 10;
    int subdivisions = 20;
    double radius = 0.075;

public:
    SpringRepresentation();
    void updateRepresentation(int index, float* positions, glm::dvec3 &a, glm::dvec3 &b);
    int getFloatCount();
};

#endif