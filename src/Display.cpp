//
//  Display.cpp
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#include "Display.h"
#include "error.h"
#include <stdio.h>


bool shaderStatus(GLuint &shader);
bool programStatus(GLuint &program);

const glm::dvec3 zaxis(0,0,1);
const glm::dvec3 xaxis(1,0,0);

int RectangularPrism::getFloatCount(){
    //N boxes with 6 sides, 2 triangles per side, 3 nodes per triangle, 3 position and 3 normal per node.
    // N          *6       *2                    *3                   * (3              +3) =
    return N*6*2*3*(3 + 3);
}

int RectangularPrism::getPositionOffset(){
    return N*SIDES*TRIANGLES*NODES*POSITIONS;
}

void RectangularPrism::updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c){
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

    glm::dvec3 norm = glm::cross(b - a, c - a);
    norm = glm::normalize(norm);
    target[0 + position_offset] = static_cast<float>(norm[0]);
    target[1 + position_offset] = static_cast<float>(norm[1]);
    target[2 + position_offset] = static_cast<float>(norm[2]);
    target[3 + position_offset] = static_cast<float>(norm[0]);
    target[4 + position_offset] = static_cast<float>(norm[1]);
    target[5 + position_offset] = static_cast<float>(norm[2]);
    target[6 + position_offset] = static_cast<float>(norm[0]);
    target[7 + position_offset] = static_cast<float>(norm[1]);
    target[8 + position_offset] = static_cast<float>(norm[2]);

}

void RectangularPrism::updateFace(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &d){
    int HALF = NODES*POSITIONS;
    updateTriangle(target, a, b, c);
    updateTriangle(&target[HALF], a, c, d);

}

int RectangularPrism::getElementNodeCount(){
    return NODES*TRIANGLES*SIDES;
}

void RectangularPrism::updateRod(int index, Rod &rod){
    int FACE = TRIANGLES*NODES*POSITIONS;

    //index times the number of floats for a rod.
    int start = index*SIDES*FACE;
    //step one create the 3 axis for changes
    glm::dvec3 main_axis(
            rod.direction[0]*0.5*(rod.length + rod.diameter),
            rod.direction[1]*0.5*(rod.length + rod.diameter),
            rod.direction[2]*0.5*(rod.length + rod.diameter)
    );

    glm::dvec3 first_axis = glm::cross(main_axis, zaxis);
    if(glm::length(first_axis)==0){
        first_axis = glm::cross(main_axis, xaxis);
    }
    first_axis = glm::normalize(first_axis);
    first_axis[0] = first_axis[0]*0.5*rod.diameter;
    first_axis[1] = first_axis[1]*0.5*rod.diameter;
    first_axis[2] = first_axis[2]*0.5*rod.diameter;

    glm::dvec3 second_axis = glm::cross(main_axis, first_axis);
    second_axis = glm::normalize(second_axis);
    second_axis[0] = second_axis[0]*0.5*rod.diameter;
    second_axis[1] = second_axis[1]*0.5*rod.diameter;
    second_axis[2] = second_axis[2]*0.5*rod.diameter;
    {
        //front
        glm::dvec3 a = rod.position + main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position + main_axis - first_axis + second_axis;
        glm::dvec3 c = rod.position + main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position + main_axis + first_axis - second_axis;

        updateFace(&positions[start], a, b, c, d);
    }

    {
        //back
        glm::dvec3 a = rod.position - main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position - main_axis + first_axis - second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position - main_axis - first_axis + second_axis;

        updateFace(&positions[start + FACE], a, b, c, d);
    }

    {
        //right
        glm::dvec3 a = rod.position + main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position + main_axis + first_axis - second_axis;
        glm::dvec3 c = rod.position - main_axis + first_axis - second_axis;
        glm::dvec3 d = rod.position - main_axis + first_axis + second_axis;

        updateFace(&positions[start + 2*FACE], a, b, c, d);
    }

    {
        //left
        glm::dvec3 a = rod.position + main_axis - first_axis + second_axis;
        glm::dvec3 b = rod.position - main_axis - first_axis + second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position + main_axis - first_axis - second_axis;

        updateFace(&positions[start + 3*FACE], a, b, c, d);
    }

    {
        //top
        glm::dvec3 a = rod.position + main_axis + first_axis + second_axis;
        glm::dvec3 b = rod.position - main_axis + first_axis + second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis + second_axis;
        glm::dvec3 d = rod.position + main_axis - first_axis + second_axis;

        updateFace(&positions[start + 4*FACE], a, b, c, d);
    }

    {
        //bottom
        glm::dvec3 a = rod.position + main_axis + first_axis - second_axis;
        glm::dvec3 b = rod.position + main_axis - first_axis - second_axis;
        glm::dvec3 c = rod.position - main_axis - first_axis - second_axis;
        glm::dvec3 d = rod.position - main_axis + first_axis - second_axis;

        updateFace(&positions[start + 5*FACE], a, b, c, d);
    }

    /*printf("box \n\n");
    for(int i = 0; i<36; i++){
        float* p = &positions[3*i];
        float* n = &positions[3*i + position_offset];
        printf("%1.1f, %1.1f, %1.1f, ... %1.1f, %1.1f, %1.1f \n", p[0], p[1], p[2], n[0], n[1], n[2]);
    }*/
}

SpringRepresentation::SpringRepresentation(){
    floats = 3*(2*rings*subdivisions + 2);
}

void SpringRepresentation::updateRepresentation(int index, float *positions, glm::dvec3 &a, glm::dvec3 &b) {
    int offset = index*floats;
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

    int n = rings*subdivisions;
    double ds = 1.0/n;
    double dtheta = 2*PI/subdivisions;
    for(int i = 0; i<2*n; i++){
        int stop = i/2;
        double theta = dtheta*stop;
        double s = ds*stop;
        positions[offset + 3*(i+1)] = (float)(a[0] + r[0]*s + a1_hat[0]*radius*cos(theta) + a2_hat[0]*radius*sin(theta));
        positions[offset + 3*(i+1) + 1] = (float)(a[1] + r[1]*s + a1_hat[1]*radius*cos(theta) + a2_hat[1]*radius*sin(theta));
        positions[offset + 3*(i+1) + 2] = (float)(a[2] + r[2]*s + a1_hat[2]*radius*cos(theta) + a2_hat[2]*radius*sin(theta));
    }

    positions[offset + 3*(2*n + 1) + 0] = (float)b[0];
    positions[offset + 3*(2*n + 1) + 1] = (float)b[1];
    positions[offset + 3*(2*n + 1) + 2] = (float)b[2];

}

int SpringRepresentation::getFloatCount(){
    return floats;
}


Display* main_display;
void keyPressedStatic(GLFWwindow* window, int key, int scancode, int action, int mods){
    main_display->keyPressed(window, key, scancode, action, mods);
};

void mousePressedStatic(GLFWwindow* window, int button, int action, int mods){
    if(action==GLFW_PRESS) {
        main_display->mousePressed(window, button, mods);
    } else{
        main_display->mouseReleased(window, button, mods);
    }
}

void mouseMovedStatic(GLFWwindow* window, double x, double y){
    main_display->mouseMoved(window, x, y);
}

Display::Display(int rods){
    //repr = new RectangularPrism(rods);
    repr = new MeshCylinder(rods, 16);
    int floats = repr->getFloatCount();

    positions = new float[floats];
    repr->setPositions(positions);

    max_springs = 100;
    current_springs = 0;
    spring_repr = new SpringRepresentation();
    spring_positions = new float[max_springs*spring_repr->getFloatCount()];


    N = rods;
    main_display=this;
}



int Display::initialize(){
    /*
     *Creating and initlizing a window on mac.
     *
     **/
    
    /* Initialize the library */
    if (!glfwInit())
        return -1;
    
#ifdef __APPLE__
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#else
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#endif

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(800, 600, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    
#ifndef __APPLE__
    glewExperimental = GL_TRUE;
    glewInit();
#endif

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);

    glEnable(GL_DEPTH_TEST);

    glDepthFunc(GL_LEQUAL);
    glDepthMask(GL_TRUE);
    glDepthRange(0.0f, 1.0f);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    printf("GLSL version %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
    printf("GL Version: %s\n", glGetString(GL_VERSION));
    
    
    /*
     * Create the program.
     */
    std::ifstream vertexFile("shaders/world.vert");
    long start = vertexFile.tellg();
    vertexFile.seekg(0, std::ios::end);
    long end = vertexFile.tellg();
    char* vertexCStr = new char[1+end-start];
    vertexCStr[end-start] = 0;
    vertexFile.seekg(0, std::ios::beg);
    vertexFile.read(vertexCStr, end-start);
    //const char* vertexCStr = vertexStr.c_str();
    GLuint vertexShader = glCreateShader( GL_VERTEX_SHADER );

    const char* vcs(vertexCStr);
    glShaderSource(vertexShader, 1, &vcs, NULL);
    glCompileShader( vertexShader );
    
    if(!shaderStatus(vertexShader)){
        exit(1);
    }


    printf("made a window\n");

    std::ifstream fragFile;
    fragFile.open("shaders/color.frag");
    start = fragFile.tellg();
    fragFile.seekg(0, std::ios::end);
    end = fragFile.tellg();
    char* fragCStr = new char[1+end-start];
    fragCStr[end-start] = 0;

    fragFile.seekg(0, std::ios::beg);
    fragFile.read(fragCStr, end-start);
    GLuint fragmentShader = glCreateShader( GL_FRAGMENT_SHADER );

    //gcc's strange complaint.
    const char* fcs(fragCStr);
    glShaderSource(fragmentShader, 1, &fcs, NULL);
    glCompileShader( fragmentShader );


    if(!shaderStatus(fragmentShader)){
        exit(1);
    }

    program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    
    if(!programStatus(program)){
        exit(1);
    }

    delete[] vertexCStr;
    vertexFile.close();
    delete[] fragCStr;
    fragFile.close();


    /*
     * Set up the buffers.
     */
    glUseProgram(program);
    
    
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    
    
    glGenBuffers(1, &positionBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, positionBufferObject);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*repr->getFloatCount(),positions, GL_STREAM_DRAW);
    
    GLuint posIndex = glGetAttribLocation(program, "position");
    GLuint normIndex = glGetAttribLocation(program, "normal");
    glEnableVertexAttribArray(posIndex);
    glVertexAttribPointer(posIndex, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(normIndex);
    glVertexAttribPointer(normIndex, 3, GL_FLOAT, GL_TRUE, 0, (GLvoid*)(repr->getPositionOffset()*sizeof(float)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);

    glGenVertexArrays(1, &vao2);
    glBindVertexArray(vao2);


    glGenBuffers(1, &springPositionBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, springPositionBufferObject);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*spring_repr->getFloatCount()*max_springs,spring_positions, GL_STREAM_DRAW);

    glEnableVertexAttribArray(posIndex);
    glVertexAttribPointer(posIndex, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);


    glUseProgram(0);
    
    GetError();

    camera = new Camera(program);

    glfwSetKeyCallback(window, keyPressedStatic);
    glfwSetMouseButtonCallback(window, mousePressedStatic);
    glfwSetCursorPosCallback(window, mouseMovedStatic);
    return 0;
}

void Display::startWriter(){
    //should be moved.
    writer = new TiffWriter("testing.tiff",height, width);
    pixbuf=new char[height*width*3];
    writing=true;
}

float myosin_color[] = {0,0,1,1};
float actin_color[] = {0.5,1,0.5,1};
float linker_color[] = {1,0,0,1};

int Display::render(){
    std::lock_guard<std::mutex> lock(mutex);
    glUseProgram(program);

        glBindBuffer(GL_ARRAY_BUFFER, positionBufferObject);
        int floats = repr->getFloatCount();
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*floats,positions);

        /* Loop until the user closes the window */
        //while
        //{
        /* Render here */
        glClearColor(0.9f, 0.9f, 1.0f, 0.0f);
        glClearDepth(1.0f);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glBindVertexArray(vao);
        GetError();

        GLint shift_loc = glGetUniformLocation(program, "shift");
        GLint color_loc = glGetUniformLocation(program, "color");
        GLint trans_loc = glGetUniformLocation(program, "transparency");
        GLint mode_loc = glGetUniformLocation(program, "colorMode");


        float* shift = new float[3];
        shift[0] = 0;shift[1]=0;shift[2]=0;shift[3]=0;
        for(int i=1; i<2; i++) {
            for(int j=1; j<2; j++) {
                shift[0] = (float)((i-1)*Constants::WIDTH);
                shift[1] = (float)((j-1)*Constants::WIDTH);
                glUniform3fv(shift_loc, 1, shift);

                if(i==1 && j==1){
                    glUniform1f(trans_loc, 1.0);
                } else{
                    glUniform1f(trans_loc, 0.25);
                 }
                glUniform1i(mode_loc, 0);
                glUniform4fv(color_loc, 1, actin_color);

                int actin_nodes = Constants::ACTINS * repr->getElementNodeCount();
                int myosin_nodes = Constants::MYOSINS * repr->getElementNodeCount();
                if (actin_nodes > 0) {
                    glDrawArrays(GL_TRIANGLES, 0, actin_nodes);
                }
                glUniform4fv(color_loc, 1, myosin_color);

                if (myosin_nodes > 0) {
                    glDrawArrays(GL_TRIANGLES, actin_nodes, myosin_nodes);
                }
            }
        }



        glBindVertexArray(0);

        if(current_springs>0) {

            glBindBuffer(GL_ARRAY_BUFFER, springPositionBufferObject);

            int spring_floats = spring_repr->getFloatCount() * max_springs;
            if(current_springs<=max_springs) {
                glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * spring_floats, spring_positions);
            } else{
                glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * spring_floats, spring_positions);
                //max_springs = current_springs;
                //glBufferData(GL_ARRAY_BUFFER, sizeof(float)*spring_repr->getFloatCount()*max_springs,spring_positions, GL_STREAM_DRAW);
            }
            glBindVertexArray(vao2);
            shift[0] = 0;
            shift[1] = 0;
            shift[2] = 0;
            shift[3] = 0;

            glUniform4fv(color_loc, 1, linker_color);
            glUniform3fv(shift_loc, 1, shift);
            glUniform1f(trans_loc, 1.0);
            glUniform1i(mode_loc, 1);

            glDrawArrays(GL_LINES, 0, spring_floats);

            glBindVertexArray(0);
        }

        glUseProgram(0);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        
        /* Poll for and process events */
        glfwPollEvents();
        if(writing) {
            glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixbuf);
            writer->writeFrame(pixbuf);
            if(writer->getCount()==last){
                writer->close();
                return -1;
            }

        }

        if(snapshot){

            takeSnapShot();

        }
    delete shift;

    return running;
    

}

void Display::shutdown(){
    if(writing) writer->close();
    glfwTerminate();
}

void Display::updateRod(int index, Rod &rod){
    //std::lock_guard<std::mutex> lock(mutex);
    repr->updateRod(index, rod);
}

void Display::setSpringCount(int s) {
    if(s>max_springs){
        std::lock_guard<std::mutex> lock(mutex);
        delete spring_positions;
        spring_positions = new float[spring_repr->getFloatCount()*s];
    }

    current_springs = s;
}

void Display::graphicsLoop(){

    while(running==0){
        render();
    }
}



int Display::getRunning(){
    return running;
}

void Display::takeSnapShot(){
    char* buf = new char[3*width*height];
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

    char* str_buf = new char[100];
    auto start = std::chrono::system_clock::now();

    sprintf(str_buf, "snapshot-%lld.tif", std::chrono::duration_cast<std::chrono::milliseconds>(start.time_since_epoch()).count());

    TiffWriter writes(str_buf, height, width);
    writes.writeFrame(buf);
    writes.close();

    delete[] str_buf;
    delete buf;
    snapshot=false;

}


void Display::keyPressed(GLFWwindow* window, int key, int scancode, int action, int mods){

    if(true){
        switch(key){
            case GLFW_KEY_LEFT:
                //camera->rotate(-0.01f, 0);
                camera->pan(0.1, 0);
                break;
            case GLFW_KEY_RIGHT:
                camera->pan(-0.1, 0);
                //camera->rotate(0.01f,0);
                break;
            case GLFW_KEY_UP:
                //camera->rotate(0, 0.01f);
                camera->pan(0, -0.1);
                break;
            case GLFW_KEY_DOWN:
                camera->pan(0, 0.1f);
                break;
            case GLFW_KEY_Z:
                camera->zoom(0.1f);
                break;
            case GLFW_KEY_A:
                camera->zoom(-0.1f);
                break;
            case GLFW_KEY_ESCAPE:
                running=-1;
                break;
            case GLFW_KEY_SPACE:
                snapshot=true;
                break;
            case GLFW_KEY_ENTER:
                releaseTrigger();
                break;
            case GLFW_KEY_G:
                moveLights(-0.1f, 0.f, 0.f);
                break;
            case GLFW_KEY_H:
                moveLights(0.1f, 0.f, 0.f);
                break;
            case GLFW_KEY_Y:
                moveLights(0.0f, 0.1f, 0.f);
                break;
            case GLFW_KEY_B:
                moveLights(0.0f, -0.1f, 0.f);
                break;
            case GLFW_KEY_N:
                moveLights(0.f, 0.f, -0.1f);
                break;
            case GLFW_KEY_J:
                moveLights(0.f, 0.f, 0.1f);

        }
    }
}



/**
 * @brief shaderStatus
 * For checking if the shader compiled and printing any error messages.
 *
 * @param shader a shader that was compiled or attempted.
 * @return true if the shader didn't fail.
 */
bool shaderStatus(GLuint &shader){
    GLint status;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE)
    {
        GLint infoLogLength;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
        
        GLchar *strInfoLog = new GLchar[infoLogLength + 1];
        glGetShaderInfoLog(shader, infoLogLength, NULL, strInfoLog);
        
        const char *strShaderType = "shader";
        
        fprintf(stderr, "Compile failure in %s shader:\n%s\n", strShaderType, strInfoLog);
        delete[] strInfoLog;
        return false;
    }
    
    
    return true;
}

/**
 * @brief programStatus
 *
 * Checks if the program linked ok.
 *
 * @param program
 * @return
 */
bool programStatus(GLuint &program){
    GLint status;
    glGetProgramiv (program, GL_LINK_STATUS, &status);
    if(status==GL_FALSE){
        GLint infoLogLength;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);
        
        GLchar *strInfoLog = new GLchar[infoLogLength + 1];
        glGetProgramInfoLog(program, infoLogLength, NULL, strInfoLog);
        fprintf(stderr, "Linker failure: %s\n", strInfoLog);
        return false;
    }
    return true;
}

MeshCylinder::MeshCylinder(int rods, int divisions){
    //        stalk points.         conical caps
    floats = rods*divisions*2*3*6 + 2*rods*divisions*3*(3+3);
    position_offset = floats/2;
    element_node_count = divisions*2*3 + 2*divisions*3;
    this->divisions = divisions;
}

int MeshCylinder::getFloatCount(){
    return floats;
}

int MeshCylinder::getPositionOffset(){
    return position_offset;
}

int MeshCylinder::getElementNodeCount() {
    return element_node_count;
}

void MeshCylinder::updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &na, glm::dvec3 &nb, glm::dvec3 &nc){
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
glm::dvec3 sum(double a, glm::dvec3 &va, double b, glm::dvec3 &vb){

    return glm::dvec3(
        a*va[0] + b*vb[0],
        a*va[1] + b*vb[1],
        a*va[2] + b*vb[2]
    );


}
void MeshCylinder::updateRod(int index, Rod &rod) {

    //index times the number of floats for a rod.
    int start = index*element_node_count*3;
    //step one create the 3 axis for changes
    glm::dvec3 main_axis(
            rod.direction[0]*0.5*(rod.length),
            rod.direction[1]*0.5*(rod.length),
            rod.direction[2]*0.5*(rod.length)
    );

    glm::dvec3 first_axis = glm::cross(main_axis, zaxis);
    if(glm::length(first_axis)==0){
        first_axis = glm::cross(main_axis, xaxis);
    }
    first_axis = glm::normalize(first_axis);

    glm::dvec3 second_axis = glm::cross(main_axis, first_axis);
    second_axis = glm::normalize(second_axis);

    double dtheta = 2*3.14159/divisions;

    double r = 0.5* rod.diameter;
    int floats_per_face = 3*6 + 6*3;

    glm::dvec3 tip = sum(0.5*rod.diameter + 0.5*rod.length, rod.direction, 1, rod.position);
    glm::dvec3 tail = sum(-0.5*(rod.diameter + rod.length), rod.direction, 1, rod.position);
    glm::dvec3 backwards(-rod.direction[0], -rod.direction[1], -rod.direction[2]);
    for(int i = 0;i<divisions; i++){
        double a1 = cos(dtheta*i);
        double b1 = sin(dtheta*i);
        double a2 = cos(dtheta*((i+1)%divisions));
        double b2 = sin(dtheta*((i+1)%divisions));

        glm::dvec3 n1 = sum(a1, first_axis, b1, second_axis);
        glm::dvec3 n2 = sum(a2, first_axis, b2, second_axis);

        glm::dvec3 a = sum(-1, main_axis, r, n1);
        glm::dvec3 b = sum(-1, main_axis, r, n2);
        glm::dvec3 c = sum(1, main_axis, r, n2);
        glm::dvec3 d = sum(1, main_axis, r, n1);

        a = sum(1, rod.position, 1, a);
        b = sum(1, rod.position, 1, b);
        c = sum(1, rod.position, 1, c);
        d = sum(1, rod.position, 1, d);

        updateTriangle(&positions[start + i*floats_per_face], a, b, d, n1, n2, n1);
        updateTriangle(&positions[start + i*floats_per_face + 9], b, c, d, n2, n2, n1);

        updateTriangle(&positions[start + i*floats_per_face + 18], d, c, tip, n1, n2, rod.direction);
        updateTriangle(&positions[start + i*floats_per_face + 27], b, a, tail, n2, n1, backwards);
    }

}


void Display::updateSpring(int index, glm::dvec3 &a, glm::dvec3 &b) {
    if(index>current_springs){
        printf("broken!\n");
    }
    spring_repr->updateRepresentation(index, spring_positions, a, b);
}

void Display::mousePressed(GLFWwindow *window, int button, int mod) {
    if(button==GLFW_MOUSE_BUTTON_LEFT){
        dragging=true;
        glfwGetCursorPos(window, &cursor_x, &cursor_y);
    }
}

void Display::mouseReleased(GLFWwindow *window, int button, int mod) {
    if(button==GLFW_MOUSE_BUTTON_LEFT){
        dragging=false;
    }
}

void Display::mouseMoved(GLFWwindow *window, double x, double y){
    if(dragging){
        double delta_x = x - cursor_x;
        double delta_y = y - cursor_y;
        camera->rotate((float)(delta_x*RATE), (float)(delta_y*RATE));
        cursor_x = x;
        cursor_y = y;
    }
}

void Display::setTrigger(std::mutex *m, std::condition_variable *cv, bool *ready) {
    starter = m;
    condition = cv;
    when_ready = ready;
}

void Display::releaseTrigger() {
    *when_ready = true;
    //std::lock_guard<std::mutex> lk(*starter);
    condition->notify_all();

}

void Display::moveLights(float dx, float dy, float dz) {
    camera->moveLight(dx, dy, dz);

}
