#ifndef _CELL_CAM
#define _CELL_CAM
#ifndef __APPLE__
    #include "glad/glad.h"
    #include "GLFW/glfw3.h"
#else
    #define GLFW_INCLUDE_GLCOREARB
    #define GLFW_NO_GLU
    #include "GLFW/glfw3.h"
#endif

class Camera {
    
    GLuint theProgram;
    float perspectiveMatrix[16];
    float orientationMatrix[16];
    float normalMatrix[9];
    
    float quarternion[4];
    float aspect;
    float x,y,z;
    float* xx;
    float* yy;
    float* zz;
    float r;
    void setPosition();
    void setRotation();
    void rotate(float* quat);
    float* light_position;
    float* ambient_light;
    float* light_intensity;
    float* location;
    public:
    
        Camera(GLuint &program);
        void resizeWindow(float w, float h);
        void rotate(float dtheta, float dphi);
        void setPerspectiveMatrix();
        void zoom(double dr);
        void pan(float dx, float dy);
        void updatePosition();
        void updateLights();
        void moveLight(float dx, float dy, float dz);
};



#endif
