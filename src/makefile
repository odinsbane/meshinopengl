CC=g++
INCLUDES=-I$(GLFW3_HOME)/include/
CFLAGS=-std=c++11 -O3
BUILD_DIR=../build
OBJDIR=$(BUILD_DIR)/obj
_OBJ = Ball.o ExecutionService.o Display.o main.o
OBJ = $(patsubst %, $(OBJDIR)/%,$(_OBJ))
DEPS = ExecutionService.h Display.h Ball.h
LIBS=-L$(GLFW3_HOME)/lib/ -lglfw3 -lrt -lXrandr -lXi -lGL -lm -ldl -lXrender -ldrm -lXdamage -lX11-xcb -lxcb-glx -lxcb-dri2 -lXxf86vm -lXfixes -lXext -lX11 -lpthread -lxcb -lXau -lGLEW
cmd_line: $(OBJ)
	$(CC) $(OBJ) $(LIBS) -o $(BUILD_DIR)/parallel_balls
	
$(OBJDIR)/%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(INCLUDES)