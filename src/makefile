CC=g++
INCLUDES=-I$(GLFW3_HOME)/include/
CFLAGS=-std=c++11 -O3
BUILD_DIR=../build
OBJDIR=$(BUILD_DIR)/obj
_OBJ = Display.o interactions.o  MeshCylinder.o MeshMyosin.o RectangularPrism.o SpringRepresentation.o Camera.o ExecutionService.o main.o MeshHelix.o rod.o Simulation.o TiffWriter.o
OBJ = $(patsubst %, $(OBJDIR)/%,$(_OBJ))
DEPS = *.h
DEPS += /usr/include/tiffio.h
LIBS=-L$(GLFW3_HOME)/lib/ -ltiff -lglfw3 -lGLEW -lrt -lm -ldl -lXrandr -lXinerama -lXi -lXcursor -lXrender -lGL -lm -lpthread -ldl -ldrm -lXdamage -lXfixes -lX11-xcb -lxcb-glx -lxcb-dri2 -lxcb-dri3 -lxcb-present -lxcb-randr -lxcb-xfixes -lxcb-render -lxcb-shape -lxcb-sync -lxshmfence -lXxf86vm -lXext -lX11 -lpthread -lxcb -lXau
cmd_line: $(OBJ)
	$(CC) $(OBJ) $(LIBS) -o $(BUILD_DIR)/parallel_balls
	
$(OBJDIR)/%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(INCLUDES)
