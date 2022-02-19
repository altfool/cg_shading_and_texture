#version 410


// Define INPUTS from fragment shader
uniform mat4 view_mat;

// And from the uniform outputs for the textures setup in main.cpp.
uniform sampler2D texture00;
uniform sampler2D texture01;

uniform vec3 light_position;
uniform float shininess;
uniform float diffuse_effect_flag;
uniform float specular_effect_flag;
uniform float use_texture_flag;
uniform float use_2nd_texture_flag;

in vec3 ex_position, ex_normal;
in vec2 ex_texture;

out vec4 fragment_color; //RGBA color

void main () {
vec3 P, N, L, V, H, ex_2nd_normal;
vec4 tcolor, tcolor2;
vec3 eye_position = vec3(0.0f, 0.0f, 5.0f);
//vec3 light_position = vec3(0.0f, 0.0f, 5.0f);
vec4 diffuse_color = vec4(0.1,0.9,0.0,1.0); 
vec4 specular_color = vec4(0.9,0.1,0.0,1.0);
vec4 ambient_color = vec4(0.1,0.1,0.1,1.0); 
//float shininess = 2.0;

tcolor2 = texture2D(texture01, ex_texture);
ex_2nd_normal = (1-use_2nd_texture_flag) * ex_normal + use_2nd_texture_flag*vec3(tcolor2[0],tcolor2[1],tcolor2[2]);
P = ex_position;
N = normalize(ex_2nd_normal);
L = normalize(light_position - P);
V = normalize(eye_position - P);
H = normalize(L+V);

tcolor = texture2D(texture00, ex_texture);
//diffuse_color = 0.1*diffuse_color + 0.9*tcolor;
diffuse_color = (1-use_texture_flag)*diffuse_color + use_texture_flag*tcolor;
diffuse_color *= max(dot(N,L),0.0);
specular_color *= pow(max(dot(H,N),0.0),shininess);
ambient_color *= diffuse_color;
//fragment_color = ambient_color;
fragment_color = ambient_color + diffuse_effect_flag*diffuse_color + specular_effect_flag*specular_color;

 
// See the normals:
//fragment_color = vec4(N,1.0); 
//fragment_color = vec4(1.0,0.5,0.0,1.0);
}
