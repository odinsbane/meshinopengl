#version 150

//input
in vec3 position;
in vec3 normal;

uniform vec3 cam_offset;
uniform mat4 perspectiveMatrix;
uniform mat4 orientationMatrix;
uniform mat3 normalModelToCameraMatrix;

uniform vec3 lightPos;
uniform vec4 lightIntensity;
uniform vec4 ambientIntensity;

smooth out float incidenceCos;
smooth out vec4 meshColor;
smooth out vec3 planePosition;

uniform vec4 color;

uniform vec3 shift;

uniform int colorMode;

void main() {

        vec4 camPosition = orientationMatrix*(vec4(position+cam_offset + shift,1));
        //dirToLight = vec3(0,0,1);

        gl_Position = perspectiveMatrix*camPosition;

        if(colorMode==0){
            vec3 dirToLight = normalize(lightPos - vec3(position));
            incidenceCos = dot(normal, dirToLight);

        }else{
            incidenceCos = -2;
        }

        meshColor = color;

}
