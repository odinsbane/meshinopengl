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
smooth out vec3 t_norm;
smooth out vec3 t_pos;

uniform vec4 color;

uniform vec3 shift;

uniform int colorMode;
uniform int toggle = 1;
void main() {
        //position in camera space.
        vec4 camPosition = orientationMatrix*(vec4(position+camOffset + shift,1));
        //dirToLight = vec3(0,0,1);

        //position in clip space.
        gl_Position = perspectiveMatrix*camPosition;
        pos = position;
        norm = normal;
        t_norm = normalize(normalModelToCameraMatrix*normal);
        t_pos = camPosition.xyz;

        if(pos.x<-3.2||pos.x>3.2||pos.y<-3.2||pos.y>3.2){
            gl_Position.z = 10*toggle + (toggle-1)*gl_Position.z;
        } else{
            gl_Position.z = 10*(1 - toggle) + toggle*gl_Position.z;
        }



        meshColor = color;


}
