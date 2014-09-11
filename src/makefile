CC=g++
INCLUDES=-I$(GLFW3_HOME)/include/
CFLAGS=-std=c++11 -O3
LIBS=-L/home/msmith/local/lib/ -lglfw3 -lrt -lXrandr -lXi -lGL -lm -ldl -lXrender -ldrm -lXdamage -lX11-xcb -lxcb-glx -lxcb-dri2 -lXxf86vm -lXfixes -lXext -lX11 -lpthread -lxcb -lXau -lGLEW
cmd_line: main.o Display.o ExecutionService.o Ball.o
	$(CC) Ball.o ExecutionService.o Display.o main.o $(LIBS) -o parallel_balls
	
Ball.o: Ball.cpp
	$(CC) $(CFLAGS) -c Ball.cpp -o Ball.o

ExecutionService.o: ExecutionService.cpp
	$(CC) $(CFLAGS) -c ExecutionService.cpp -o ExecutionService.o

Display.o: Display.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c Display.cpp -o Display.o

main.o: main.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c main.cpp -o main.o
