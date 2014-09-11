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
    char* pixbuf;
    int height = 300;
    int width = 400;
public:
    Display(int N);
    int initialize();
    void updateBall(int index, double x, double y, double radius);
    int render();
    ~Display(){
        
    }
};
#endif /* defined(__ParallelBalls__Display__) */
