#!/bin/bash
#script for building using pkg-config
#the argmument passed should be the glfw3 package config file.
export GLFW_PC=$1
export GLFW_CFLAGS=$(pkg-config --cflags "$GLFW_PC")
export GLFW_LIBS=$(pkg-config --static --libs "$GLFW_PC")
export GLAD_HOME=$2

g++ -std=c++11 -I$GLAD_HOME/include $GLFW_CFLAGS src/*.cpp $GLAD_HOME/src/glad.c $GLFW_LIBS -lGL -lglapi -lGLdispatch -ltiff -o build/fastshadows

