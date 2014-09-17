//
//  Display.cpp
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#include "Display.h"
#include "error.h"



bool shaderStatus(GLuint &shader);
bool programStatus(GLuint &program);

Display::Display(int nodes){
    //N boxes with 6 sides, 2 triangles per side, 3 nodes per triangle, 3 position and 3 normal per node.
    // N          *6       *2                    *3                   * (3              +3) =
    positions = new float[6*2*3*(3+3)*nodes];
    N = nodes;
    position_offset = 6*2*3*3;
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
#endif
    
    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(400, 300, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    
#ifndef __APPLE__
        glewInit();
#endif
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
    char* vertexCStr = new char[end-start];
    vertexFile.seekg(0, std::ios::beg);
    vertexFile.read(vertexCStr, end-start);
    //const char* vertexCStr = vertexStr.c_str();
    GLuint vertexShader = glCreateShader( GL_VERTEX_SHADER );
    glShaderSource(vertexShader, 1, &vertexCStr, NULL);
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
    char* fragCStr = new char[end-start];

    fragFile.seekg(0, std::ios::beg);
    fragFile.read(fragCStr, end-start);
    GLuint fragmentShader = glCreateShader( GL_FRAGMENT_SHADER );
    glShaderSource(fragmentShader, 1, &fragCStr, NULL);
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
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*72*N,positions, GL_STREAM_DRAW);
    
    GLuint posIndex = glGetAttribLocation(program, "position");
    GLuint normIndex = glGetAttribLocation(program, "normal");
    glEnableVertexAttribArray(posIndex);
    glVertexAttribPointer(posIndex, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(normIndex);
    glVertexAttribPointer(normIndex, 3, GL_FLOAT, GL_TRUE, position_offset, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);
    glUseProgram(0);
    
    GetError();

    //should be moved.
    writer = new TiffWriter("testing.tiff",height, width);
    pixbuf=new char[height*width*3];

    camera = new Camera(program);

    return 0;
}

int Display::render(){
        glUseProgram(program);
        glBindBuffer(GL_ARRAY_BUFFER, positionBufferObject);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*4*6*N,positions, GL_STREAM_DRAW);
    
        /* Loop until the user closes the window */
        //while
        //{
        /* Render here */
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        glClearDepth(1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUseProgram(program);
        glBindVertexArray(vao);
        
        glDrawArrays(GL_TRIANGLES, 0, 6*N);
        
        glBindVertexArray(0);
        glUseProgram(0);
        
        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        
        /* Poll for and process events */
        glfwPollEvents();
        if(writer->isOpen()) {
            glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixbuf);
            writer->writeFrame(pixbuf);
            if(writer->getCount()==last){
                writer->close();
                return -1;
            }

        }
        if(glfwWindowShouldClose(window)) return -1;
    
        return 0;
    

}

void Display::shutdown(){
    writer->close();
    glfwTerminate();
}

const glm::dvec3 zaxis(0,0,1);
const glm::dvec3 xaxis(1,0,0);

void Display::updateRod(int index, Rod &rod){
    int start = index*72;
    //step one create the 3 axis for changes
    glm::dvec3 main_axis(
            rod.direction[0]*0.5*rod.length,
            rod.direction[1]*0.5*rod.length,
            rod.direction[2]*0.5*rod.length
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

    glm::dvec3 a = rod.position + main_axis + first_axis + second_axis;
    glm::dvec3 b = rod.position + main_axis + first_axis - second_axis;
    glm::dvec3 c = rod.position + main_axis - first_axis - second_axis;
    glm::dvec3 d = rod.position + main_axis - first_axis + second_axis;

    updateFace(&positions[start], a, b, c, d);

}

void Display::updateFace(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c, glm::dvec3 &d){

    updateTriangle(target, a, b, c);
    updateTriangle(target + 36, b, d, c);

}

void Display::updateTriangle(float* target, glm::dvec3 &a, glm::dvec3 &b, glm::dvec3 &c){
    target[0] = a[0];
    target[1] = a[1];
    target[2] = a[2];
    target[3] = b[0];
    target[4] = b[1];
    target[5] = b[2];
    target[6] = c[0];
    target[7] = c[1];
    target[8] = c[2];


    glm::dvec3 norm = glm::cross(b - a, c - a);

    target[0 + position_offset] = norm[0];
    target[1 + position_offset] = norm[1];
    target[2 + position_offset] = norm[2];
    target[3 + position_offset] = norm[0];
    target[4 + position_offset] = norm[1];
    target[5 + position_offset] = norm[2];
    target[6 + position_offset] = norm[0];
    target[7 + position_offset] = norm[1];
    target[8 + position_offset] = norm[2];
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
