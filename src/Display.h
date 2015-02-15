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
#include <condition_variable>
#include "Representations.h"

class Display{
private:
    GLFWwindow* window;

    float* positions;

    float* spring_positions;
    int max_springs;
    int current_springs;
    int a_count, m_count;


    GLuint program;
    GLuint vao;
    GLuint vao2;
    GLuint positionBufferObject;
    GLuint springPositionBufferObject;

    bool dragging = false;
    double cursor_x, cursor_y;

    TiffWriter* writer;
    char* pixbuf;
    int height = 600;
    int width = 800;
    int last = 2000;
    bool writing=false;
    bool waiting_to_write=false;
    Camera* camera;
    int running = 0;
    CylinderRepresentation *actin_repr, *myosin_repr;
    SpringRepresentation* spring_repr;
    bool snapshot=false;
    std::mutex mutex;
    std::mutex* starter;
    std::condition_variable* condition;
    bool* when_ready;
public:
    Display();
    void setRodCounts(int actins, int myosins);
    int initialize();
    void updateRod(int index, Rod &rod);
    void updateSpring(int index, glm::dvec3 &a, glm::dvec3 &b);
    void setSpringCount(int s);
    int render();
    void shutdown();
    void startWriter();
    void requestNextFrame();
    void takeSnapShot();
    void graphicsLoop();
    int getRunning();
    ~Display(){
        delete[] pixbuf;
    }
    void keyPressed(GLFWwindow* window, int key, int scancode, int action, int mods);
    void mousePressed(GLFWwindow* window, int button, int mod);
    void mouseReleased(GLFWwindow* window, int button, int mod);
    void mouseMoved(GLFWwindow* window, double x, double y);
    void updateLights();
    void setTrigger(std::mutex* m, std::condition_variable* cv, bool* ready);
    void releaseTrigger();
    double RATE = 0.001;

    void moveLights(float dx, float dy, float dz);
};


#endif /* defined(__ParallelBalls__Display__) */
