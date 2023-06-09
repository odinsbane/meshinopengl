#version 150
//version defined in program

//position ordinates are located in gl_FragCoord
//and are in 'window coordinates'
out vec4 outputColor;

uniform mat3 normalModelToCameraMatrix;

smooth in vec4 meshColor;
smooth in vec3 pos;
smooth in vec3 norm;

smooth in vec3 t_norm;
smooth in vec3 t_pos;

uniform vec4 lightIntensity;
uniform vec4 ambientIntensity;

uniform vec3 lightPos;

float transparency;
vec4 WHITE=vec4(1,1,1,1);
float SPEC = 20;
void main() {
    vec4 oc;
    
	vec3 disp = lightPos - pos;
	float l = dot(disp, disp);
	vec3 lDir = normalize(normalModelToCameraMatrix*disp);

	float incidenceCos = dot(t_norm, lDir);
	if(incidenceCos<0){
		incidenceCos=0;
	}
	vec3 n_pos = normalize(t_pos);
	vec3 rDir = reflect(-lDir, n_pos);
	vec3 m_norm = vec3(-t_norm.x, -t_norm.y, t_norm.z);
	float phong = dot(rDir, m_norm);
	phong = clamp(phong, 0, 1);
	float f = 1/pow(l/5.0, 0.25);
	
	//phong = incidenceCos!=0?phong:0;
        phong = pow(phong, SPEC)*0;



	oc = meshColor*(incidenceCos*lightIntensity*f + meshColor*ambientIntensity) + phong*WHITE;
    if(pos.x<=-3.2||pos.x>=3.2||pos.y<=-3.2||pos.y>=3.2){

        if(t_norm.y<0&&t_norm.z<0.1){
            outputColor = vec4(0,0,0, 1.0);
        } else{
            outputColor = vec4(oc.xyz, 0.5);
        }
    } else{
       outputColor = vec4(oc.xyz, 1.0);
    }



}


