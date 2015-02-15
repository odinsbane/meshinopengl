#version 150
//version defined in program

//position ordinates are located in gl_FragCoord
//and are in 'window coordinates'
out vec4 outputColor;

smooth in vec4 meshColor;
smooth in vec3 pos;
smooth in vec3 norm;
smooth in vec3 tnorm;

in float mode;

uniform float transparency;
uniform vec4 lightIntensity;
uniform vec4 ambientIntensity;

uniform vec3 lightPos;

vec4 WHITE=vec4(1,1,1,1);
float SPEC = 5;
void main() {
    if(mode==0){
        vec3 disp = lightPos - pos;
        float l = dot(disp, disp);
        vec3 lDir = normalize(disp);
        float incidenceCos = dot(norm, lDir);
        if(incidenceCos<0){
            incidenceCos=0;
        }
        vec3 cs_pos = normalize(pos);
        vec3 rDir = reflect(-lDir, cs_pos);
        float phong = dot(rDir, cs_pos);
        phong = clamp(phong, 0, 1);
        float f = 1/pow(l, 0.25);
        phong = incidenceCos!=0?phong:0;
        phong = pow(phong, SPEC);
        vec4 oc = meshColor*(incidenceCos*lightIntensity*f + meshColor*ambientIntensity);
        outputColor = vec4(oc.xyz, transparency);

    }
    if(mode==1){
        outputColor = vec4(meshColor.xyz, 1);
    }

    if(mode>1){
        outputColor=vec4(1,1,1,1);
    }

}


