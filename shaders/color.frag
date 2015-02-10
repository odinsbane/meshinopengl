#version 150
//version defined in program

//position ordinates are located in gl_FragCoord
//and are in 'window coordinates'
out vec4 outputColor;

smooth in vec4 meshColor;
smooth in vec3 pos;
smooth in vec3 norm;

in float mode;

uniform float transparency;
uniform vec4 lightIntensity;
uniform vec4 ambientIntensity;

uniform vec3 lightPos;


void main() {
    if(mode==0){
        vec3 disp = lightPos - pos;
        float l = dot(disp, disp);
        float incidenceCos = dot(norm, normalize(disp));
        vec4 oc = meshColor*(incidenceCos*lightIntensity + meshColor*ambientIntensity);
        outputColor = vec4(oc.xyz, transparency);

    }
    if(mode==1){
        outputColor = vec4(meshColor.xyz, 1);
    }

    if(mode>1){
        outputColor=vec4(1,1,1,1);
    }

}


