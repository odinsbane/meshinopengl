#version 150
//version defined in program

//position ordinates are located in gl_FragCoord
//and are in 'window coordinates'
out vec4 outputColor;

smooth in float incidenceCos;
smooth in vec4 meshColor;
smooth in vec3 planePosition;

uniform float transparency;
uniform vec4 lightIntensity;
uniform vec4 ambientIntensity;

void main() {

    if(incidenceCos==-2){
         outputColor = vec4(meshColor.xyz, transparency);
    } else if(incidenceCos>0.95){
        outputColor = vec4(1,1,1,transparency);
    }else if(incidenceCos>0){

        vec4 oc = meshColor*(incidenceCos*lightIntensity + meshColor*ambientIntensity);
        outputColor = vec4(oc.xyz, transparency);
    } else{
        vec4 oc =  meshColor*ambientIntensity;
        outputColor = vec4(oc.xyz, transparency);
    }



}


