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

smooth out vec4 meshColor;
smooth out vec3 planePosition;

vec4 color = vec4(1,0,0,1);
void main() {
        vec4 camPosition = orientationMatrix*(vec4(position+cam_offset,1));
        //dirToLight = vec3(0,0,1);

        gl_Position = perspectiveMatrix*camPosition;
        //gl_Position.z = camPosition.z;
        planePosition = vec3(0,0,0);

        vec3 nn = normalize(normal);
        vec3 dirToLight = normalize(lightPos - vec3(position));
        float cosAngIncidence = dot(nn, dirToLight)*0.5 + 0.5;
        if(color.x<color.y){
            cosAngIncidence = 1;
        }
        meshColor = (color* lightIntensity*cosAngIncidence ) + color * ambientIntensity;

}
