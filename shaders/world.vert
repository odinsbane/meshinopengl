#version 150

//input
in vec3 position;
in vec3 normal;

uniform vec3 camOffset;
uniform mat4 perspectiveMatrix;
uniform mat4 orientationMatrix;
uniform mat3 normalModelToCameraMatrix;


smooth out vec4 meshColor;
smooth out vec3 norm;
smooth out vec3 pos;
smooth out vec3 tnorm;

out float mode;

uniform vec4 color;

uniform vec3 shift;

uniform int colorMode;

void main() {

        vec4 camPosition = orientationMatrix*(vec4(position+camOffset + shift,1));
        //dirToLight = vec3(0,0,1);

        gl_Position = perspectiveMatrix*camPosition;

        if(colorMode==0){
            mode = 0;
            //norm = normalize(normalModelToCameraMatrix*normal);
            pos = position;
            norm = normal;
            tnorm = normalize(normalModelToCameraMatrix*normal);

        }else{
            mode=1;
            norm = normal;
            pos = position;
        }

        meshColor = color;

}
