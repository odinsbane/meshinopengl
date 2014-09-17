#version 150
//version defined in program

//position ordinates are located in gl_FragCoord
//and are in 'window coordinates'
out vec4 outputColor;

smooth in vec4 meshColor;
smooth in vec3 planePosition;
void main() {
    outputColor = meshColor;
}


