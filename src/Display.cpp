//
//  Display.cpp
//  ParallelBalls
//
//  Created by msmith on 9/8/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#include "Display.h"
#include "error.h"


//simple shader output = input. texture coordinates are the xy, should make 4 copies of texture.
const std::string vertexStr(
#ifdef __APPLE__
                            "#version 150\n"
#else
                            "#version 130\n"
#endif
                            "in vec4 pos;\n"
                            "out vec2 texCoords;\n"
                            "void main(){\n"
                            "   texCoords=vec2(pos.x-pos.z, pos.y-pos.w);\n"
                            "   gl_Position = vec4(pos.x, pos.y, 0, 1.0);\n"
                            "}\n"
                            );

//Gets the values from the texture.
const std::string fragmentStr(
#ifdef __APPLE__
                              "#version 150\n"
#else
                              "#version 130\n"
#endif
                              "out vec4 outputColor;\n"
                              "in vec2 texCoords;\n"
                              "uniform float limit = 3.6e-5;\n"
                              "void main(){\n"
                              "float m = dot(texCoords, texCoords);\n"
                              "if(m<limit){\n"
                              "outputColor = vec4(0,0,1,1);\n"
                              "}\n"
                              "else outputColor = vec4(0,0,0,0);\n"
                              "}\n"
                              );
bool shaderStatus(GLuint &shader);
bool programStatus(GLuint &program);

Display::Display(int nodes){
    //two positions per 6 vertices for each node
    positions = new float[4*6*nodes];
    N = nodes;

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
    const char* vertexCStr = vertexStr.c_str();
    GLuint vertexShader = glCreateShader( GL_VERTEX_SHADER );
    glShaderSource(vertexShader, 1, &vertexCStr, NULL);
    glCompileShader( vertexShader );
    
    if(!shaderStatus(vertexShader)){
        exit(1);
    }
    printf("made a window\n");
    const char* fragmentCStr = fragmentStr.c_str();
    GLuint fragmentShader = glCreateShader( GL_FRAGMENT_SHADER );
    glShaderSource(fragmentShader, 1, &fragmentCStr, NULL);
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
    
    /*
     * Set up the buffers.
     */
    glUseProgram(program);
    
    
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    
    
    glGenBuffers(1, &positionBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, positionBufferObject);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*4*6*N,positions, GL_STREAM_DRAW);
    
    GLuint posIndex = glGetAttribLocation(program, "pos");
    glEnableVertexAttribArray(posIndex);
    glVertexAttribPointer(posIndex, 4, GL_FLOAT, GL_FALSE, 0, 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);
    glUseProgram(0);
    
    GetError();

    //should be moved.
    writer = new TiffWriter("testing.tiff",height, width);
    pixbuf=new char[height*width*3];

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
        int last = 2000;
        if(writer->getCount()<last) {
            glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixbuf);
            writer->writeFrame(pixbuf);
            if(writer->getCount()==last){
                writer->close();
                std::cout<<"finished writing\n";
                glfwTerminate();
            }
        }
        if(glfwWindowShouldClose(window)) return -1;
    
        return 0;
    

}

void Display::updateBall(int index, double x, double y, double radius){
    float* node = &positions[24*index];
    node[0] = x - radius;
    node[1] = y - radius;
    node[2] = x;
    node[3] = y;
    node[4] = x + radius;
    node[5] = y - radius;
    node[6] = x;
    node[7] = y;
    node[8] = x - radius;
    node[9] = y + radius;
    node[10] = x;
    node[11] = y;
    node[12] = x + radius;
    node[13] = y - radius;
    node[14] = x;
    node[15] = y;
    node[16] = x - radius;
    node[17] = y + radius;
    node[18] = x;
    node[19] = y;
    node[20] = x + radius;
    node[21] = y + radius;
    node[22] = x;
    node[23] = y;
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
