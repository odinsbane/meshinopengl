//
//  Display.h
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//
#ifndef __APPLE__
  #include <GL/glew.h>
#else
  #define GLFW_INCLUDE_GLCOREARB
  #define GLFW_NO_GLU
#endif

#include "GLFW/glfw3.h"
#include "TiffWriter.h"
#include "rod.h"
#include "Camera.h"
#include <fstream>
#include "Simulation.h"
#ifndef __ParallelBalls__Display__
#define __ParallelBalls__Display__

#include <iostream>
#include <ctime>
#include <string>
#include <thread>



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

class SpringRepresentation{
    int divisions, floats;

public:
    SpringRepresentation();
    void updateRepresentation(int index, float* positions, glm::dvec3 &a, glm::dvec3 &b);
    int getFloatCount();
};

class Display{
private:
    GLFWwindow* window;

    float* positions;

    float* spring_positions;
    int max_springs;
    int current_springs;
    int to_set;

    int N;

    GLuint program;
    GLuint vao;
    GLuint vao2;
    GLuint positionBufferObject;
    GLuint springPositionBufferObject;

    TiffWriter* writer;
    char* pixbuf;
    int height = 300;
    int width = 400;
    int last = 2000;
    bool writing=false;
    Camera* camera;
    int running = 0;
    CylinderRepresentation* repr;
    SpringRepresentation* spring_repr;
    bool snapshot=false;
    std::mutex mutex;
public:
    Display(int N);
    int initialize();
    void updateRod(int index, Rod &rod);
    void updateSpring(int index, glm::dvec3 &a, glm::dvec3 &b);
    void setSpringCount(int s);
    int render();
    void shutdown();
    void startWriter();
    void takeSnapShot();
    void graphicsLoop();
    int getRunning();
    ~Display(){
        delete[] pixbuf;
    }
    void keyPressed(GLFWwindow* window, int key, int scancode, int action, int mods);
    void updateLights();

};


#endif /* defined(__ParallelBalls__Display__) */
