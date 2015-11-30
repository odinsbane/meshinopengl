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

float myosin_color[] = {0.0,0.0,1.0,1};
float actin_color[] = {1.0,0.0,0.0,1};
float linker_color[] = {0.33,1.0,0.33,1};
float bg[] = {1,1,1};

bool shaderStatus(GLuint &shader);
bool programStatus(GLuint &program);


Display* main_display;
void keyPressedStatic(GLFWwindow* window, int key, int scancode, int action, int mods){
    main_display->keyPressed(window, key, scancode, action, mods);
};

void windowFocusCallback(GLFWwindow* window, int focused)
{
    if (focused)
    {
        printf("gained\n");
        printf("%d ...\n",glfwGetKey(window, GLFW_KEY_A));

    }
    else
    {
        printf("lost\n");
        printf("%d ...\n",glfwGetKey(window, GLFW_KEY_A));

    }
}



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

Display::Display(){


    max_springs = 100;
    current_springs = 0;
    spring_repr = new SpringRepresentation();
    spring_repr->setMaxSpringCount(max_springs);
    spring_positions = new float[max_springs*spring_repr->getFloatCount()];


    main_display=this;
}

void Display::setRodCounts(int actins, int myosins){
    actin_repr = new MeshHelix(actins);
    int actin_floats = actin_repr->getFloatCount();
    myosin_repr = new MeshMyosin(myosins);
    int myosin_floats = myosin_repr->getFloatCount();
    int total = actin_floats + myosin_floats;
    positions = new float[total];
    actin_repr->setPositions(positions);
    myosin_repr->setPositions(positions);

    actin_repr->setPositionOffset(total/2);
    myosin_repr->setPositionOffset(total/2);
    a_count = actins;
    m_count = myosins;
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
    window = glfwCreateWindow(width, height, "Hello World", NULL, NULL);
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
    glEnable(GL_BLEND);
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

    int total_floats = actin_repr->getFloatCount() + myosin_repr->getFloatCount();
    
    glGenBuffers(1, &positionBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, positionBufferObject);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_floats,positions, GL_STREAM_DRAW);
    
    GLuint posIndex = glGetAttribLocation(program, "position");
    GLuint normIndex = glGetAttribLocation(program, "normal");

    glEnableVertexAttribArray(posIndex);
    glVertexAttribPointer(posIndex, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(normIndex);
    glVertexAttribPointer(normIndex, 3, GL_FLOAT, GL_TRUE, 0, (GLvoid*)(total_floats/2*sizeof(float)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);

    glGenVertexArrays(1, &vao2);
    glBindVertexArray(vao2);


    glGenBuffers(1, &springPositionBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, springPositionBufferObject);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*spring_repr->getFloatCount()*max_springs,spring_positions, GL_STREAM_DRAW);

    glEnableVertexAttribArray(posIndex);
    glVertexAttribPointer(posIndex, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(normIndex);
    glVertexAttribPointer(normIndex, 3, GL_FLOAT, GL_TRUE, 0, (GLvoid*)(max_springs*spring_repr->getFloatCount()/2*sizeof(float)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);


    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);


    glUseProgram(0);
    
    GetError();

    camera = new Camera(program);

    glfwSetKeyCallback(window, keyPressedStatic);
    glfwSetMouseButtonCallback(window, mousePressedStatic);
    glfwSetCursorPosCallback(window, mouseMovedStatic);
    glfwSetWindowFocusCallback(window, windowFocusCallback);

    return 0;
}

void Display::startWriter(){
    printf("recording\n");
    //should be moved.
    glfwGetFramebufferSize(window, &width, &height );
    writer = new TiffWriter("testing.tiff",height, width);
    pixbuf=new char[height*width*3];
    writing=true;
}

void Display::requestNextFrame() {
    if(!writing){
        return;
    }
    waiting_to_write=true;
}



int Display::render(){
    std::lock_guard<std::mutex> lock(mutex);
    glUseProgram(program);

        glBindBuffer(GL_ARRAY_BUFFER, positionBufferObject);
        int floats = actin_repr->getFloatCount() + myosin_repr->getFloatCount();
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*floats,positions);

        /* Loop until the user closes the window */
        //while
        //{
        /* Render here */
        glClearColor(bg[0], bg[1], bg[2], 0.0f);
        glClearDepth(1.0f);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glBindVertexArray(vao);
        GetError();

        GLint shift_loc = glGetUniformLocation(program, "shift");
        GLint color_loc = glGetUniformLocation(program, "color");
        GLint trans_loc = glGetUniformLocation(program, "transparency");
        GLint mode_loc = glGetUniformLocation(program, "colorMode");
        GLint tog_loc = glGetUniformLocation(program, "toggle");

        float* shift = new float[3];
        shift[0] = 0;shift[1]=0;shift[2]=0;shift[3]=0;
        glEnable(GL_CULL_FACE);
        for(int i=0; i<2; i++) {
            for(int j=1; j<2; j++) {
                //shift[0] = (float)((i-1)*Constants::WIDTH);
                //shift[1] = (float)((j-1)*Constants::WIDTH);
                glUniform3fv(shift_loc, 1, shift);
                glUniform1i(tog_loc, (i+1)%2);
                if(i==1 && j==1){
                    glUniform1f(trans_loc, 1.0);
                } else{
                    glUniform1f(trans_loc, 0.25);
                 }
                glUniform1i(mode_loc, 0);
                glUniform4fv(color_loc, 1, actin_color);

                int actin_nodes = a_count * actin_repr->getElementNodeCount();
                int myosin_nodes = m_count * myosin_repr->getElementNodeCount();
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
            glDisable(GL_CULL_FACE);
            glBindVertexArray(vao2);
            glBindBuffer(GL_ARRAY_BUFFER, springPositionBufferObject);


            int spring_floats = spring_repr->getFloatCount() * max_springs;
            if(current_springs<=max_springs) {
                glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * spring_floats, spring_positions);
            } else{
                //glBufferData(GL_ARRAY_BUFFER, 0, sizeof(float) * spring_floats, spring_positions);

                max_springs = current_springs;
                glBufferData(GL_ARRAY_BUFFER, sizeof(float)*spring_repr->getFloatCount()*max_springs,spring_positions, GL_STREAM_DRAW);

                GLuint normIndex = glGetAttribLocation(program, "normal");

                glEnableVertexAttribArray(normIndex);
                glVertexAttribPointer(normIndex, 3, GL_FLOAT, GL_TRUE, 0, (GLvoid*)(max_springs*spring_repr->getFloatCount()/2*sizeof(float)));
                glBindBuffer(GL_ARRAY_BUFFER, 0);


                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glBindVertexArray(0);


            }
            //glBindVertexArray(vao2);
            shift[0] = 0;
            shift[1] = 0;
            shift[2] = 0;
            shift[3] = 0;

            glUniform4fv(color_loc, 1, linker_color);
            glUniform3fv(shift_loc, 1, shift);
            glUniform1f(trans_loc, 1.0);
            glUniform1i(mode_loc, 1);

            int chunk = spring_repr->getFloatCount()/6;
            for(int i = 0; i<current_springs; i++){
                glDrawArrays(GL_TRIANGLE_STRIP, i*chunk, chunk);
            }


            glBindVertexArray(0);
        }

        glUseProgram(0);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        
        /* Poll for and process events */
        glfwPollEvents();
        if(waiting_to_write) {
            glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixbuf);
            writer->writeFrame(pixbuf);
            if(writer->getCount()==last){
                writer->close();
                return -1;
            }
            waiting_to_write=false;

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
    if(index>=a_count){

        int start = (index-a_count)*myosin_repr->getElementNodeCount()*3 + a_count*actin_repr->getElementNodeCount()*3;

        myosin_repr->updateRod(start, rod);
    }else{
        int start = index*actin_repr->getElementNodeCount()*3;
        actin_repr->updateRod(start, rod);
    }
}

void Display::setSpringCount(int s) {
    if(s>max_springs){
        std::lock_guard<std::mutex> lock(mutex);
        delete spring_positions;
        spring_positions = new float[spring_repr->getFloatCount()*s];
        spring_repr->setMaxSpringCount(s);

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
    int w,h;
    glfwGetFramebufferSize(window, &w, &h);

    char* buf = new char[3*w*h];
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, buf);

    char* str_buf = new char[100];
    auto start = std::chrono::system_clock::now();

    sprintf(str_buf, "snapshot-%lld.tif", std::chrono::duration_cast<std::chrono::milliseconds>(start.time_since_epoch()).count());

    TiffWriter writes(str_buf, h, w);
    writes.writeFrame(buf);
    writes.close();

    delete[] str_buf;
    delete buf;
    snapshot=false;

}


void Display::keyPressed(GLFWwindow* window, int key, int scancode, int action, int mods){
    printf("%d, %d, %d, %d \n", key, scancode, action, mods);
    if(action!=0){
        switch(key){
            case GLFW_KEY_LEFT:
                //camera->rotate(-0.01f, 0);
                camera->pan(0.01, 0);
                break;
            case GLFW_KEY_RIGHT:
                camera->pan(-0.01, 0);
                //camera->rotate(0.01f,0);
                break;
            case GLFW_KEY_UP:
                //camera->rotate(0, 0.01f);
                camera->pan(0, -0.01);
                break;
            case GLFW_KEY_DOWN:
                camera->pan(0, 0.01f);
                break;
            case GLFW_KEY_Z:
                printf("z pressed\n");
                camera->zoom(0.05f);
                break;
            case GLFW_KEY_A:
                printf("a pressed\n");
                //exit(-1);
                break;
            case GLFW_KEY_Q:
                printf("%d ...\n",glfwGetKey(window, GLFW_KEY_A));
                camera->zoom(-0.05f);
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
                break;
            case GLFW_KEY_R:
                startWriter();
                break;

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
