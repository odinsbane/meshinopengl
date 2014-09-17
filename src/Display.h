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

#ifndef __ParallelBalls__Display__
#define __ParallelBalls__Display__

#include <iostream>
class Display{
private:
    GLFWwindow* window;
    float* positions;
    int N;
    GLuint program;
    GLuint vao;
    GLuint positionBufferObject;
    TiffWriter* writer;
    void updateFace(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &d);
    void updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c);
    char* pixbuf;
    int height = 300;
    int width = 400;
    int last = 2000;
    int position_offset;
    Camera* camera;
public:
    Display(int N);
    int initialize();
    void updateRod(int index, Rod &rod);
    int render();
    void shutdown();
    ~Display(){
        delete[] pixbuf;
    }

    void updateLights();
};
#endif /* defined(__ParallelBalls__Display__) */
