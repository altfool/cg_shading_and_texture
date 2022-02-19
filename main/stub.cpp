/******************************************************************************|
| CPSC 4050/6050 Computer Garphics Assignment 5, Daljit Singh Dhillon, 2020    |
| Reference:                                                                   |
|                                                                              |
| Some OpenGL setup code here including math_funcs, and gl_utils               |
| are from Angton Gerdelan and "Anton's OpenGL 4 Tutorials."                   |
| http://antongerdelan.net/opengl/                                             |
| Email: anton at antongerdelan dot net                                        |
| Copyright Dr Anton Gerdelan, Trinity College Dublin, Ireland.                |
|******************************************************************************/
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <assert.h>


#include <math.h>
#include <time.h>

#include "maths_funcs.h"   // Anton's maths functions.
#include "gl_utils.h"      // Anton's opengl functions and small utilities like logs
#include "stb_image.h"     // Sean Barrett's image loader with Anton's load_texture()

#define _USE_MATH_DEFINES
#define ONE_DEG_IN_RAD (2.0 * M_PI) / 360.0 // 0.017444444
#define MY_PI 3.1415926

mat4 view_mat;
mat4 proj_mat;
mat4 model_mat;


int pointCount=0;

//forwarding definition
void binomialCoeffs(int n, int* C);
GLfloat* getRevolutionVertices(float* cpX, float* cpY, int nCtrlPts, int ysteps, int thetasteps);
GLfloat* getRevolutionNormals(float* cpX, float* cpY, int nCtrlPts, int ysteps, int thetasteps);
GLfloat* getRevolutionTexCoords(int ysteps, int thetasteps);

void loadSurfaceOfRevolution() 
{
/*------------------------------CREATE GEOMETRY-------------------------------*/
/*	GLfloat vp[18];    // array of vertex points
	
	//face 1, vertex 1
	vp[0] = -1; //x
	vp[1] = -1; //y
	vp[2] = 0; //z
	//face 1, vertex 2
	vp[3] = 1; //x
	vp[4] = -1; //y
	vp[5] = 0; //z
	//face 1, vertex 3
	vp[6] = -1; //x
	vp[7] =  1; //y
	vp[8] =  0; //z
	
	//face 2, vertex 1
	vp[ 9] = -1; //x
	vp[10] =  1; //y
	vp[11] = 0; //z
	//face 2, vertex 2
	vp[12] =  1; //x
	vp[13] = -1; //y
	vp[14] = 0; //z
	//face 2, vertex 3
	vp[15] =  1; //x
	vp[16] =  1; //y
	vp[17] =  0; //z
*/
	// VAO -- vertex attribute objects bundle the various things associated with vertices
	GLuint vao;
	glGenVertexArrays (1, &vao);   // generating and binding is common pattern in OpenGL
	glBindVertexArray (vao);       // basically setting up memory and associating it

	// VBO -- vertex buffer object to contain coordinates
	// MODIFY THE FOLLOWING BLOCK OF CODE APPRORIATELY FOR YOUR SURFACE OF REVOLUTION
	float cpX[6]={-0.2,-0.9,-0.1,-0.1,-0.9,-0.2}; float cpY[6]={-0.9,-0.7,-0.1,0.1,0.7,0.9}; 
	int nCtrlPts=6, ysteps=500, thetasteps=500;
	pointCount = 18*ysteps*thetasteps;
	// ==== task A ======
	GLfloat* myVp = getRevolutionVertices(cpX, cpY, nCtrlPts, ysteps, thetasteps);
	GLuint points_vbo;
	glGenBuffers(1, &points_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	glBufferData(GL_ARRAY_BUFFER, pointCount * sizeof (GLfloat), myVp, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);
	delete [] myVp;
	// VBO -- normals -- needed for shading calculations
	// ADD CODE TO POPULATE AND LOAD PER-VERTEX SURFACE NORMALS  
	// [HINT] Vertex normals are organized in same order as that for vertex coordinates
	// ===== task B ======
	GLfloat* myVpNormals = getRevolutionNormals(cpX, cpY, nCtrlPts, ysteps, thetasteps);
	GLuint normals_vbo;
	glGenBuffers(1, &normals_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
	glBufferData(GL_ARRAY_BUFFER, pointCount * sizeof (GLfloat), myVpNormals, GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(1);
	delete [] myVpNormals;
    	// VBO -- vt -- texture coordinates
	// ADD CODE TO POPULATE AND LOAD PER-VERTEX TEXTURE COORDINATES  
	// [HINT] texture coordinates are organized in same order as that for vertex coordinates
	// [HINT] there are two texture coordinates instead of three vertex coordinates for each vertex
	// ===== task C ======
	GLfloat* myVpTexCoords = getRevolutionTexCoords(ysteps, thetasteps);
	GLuint texture_vbo;
	glGenBuffers(1, &texture_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, texture_vbo);
	glBufferData(GL_ARRAY_BUFFER, pointCount/3*2 * sizeof (GLfloat), myVpTexCoords, GL_STATIC_DRAW);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(2);
	delete [] myVpTexCoords;
}


float light_x=0, light_y=0, light_z=5;
float light_change=0.2;
float shininess_value=2, shininess_change=0.1;
float diffuse_effect_flag=1;
float specular_effect_flag=1;
float use_texture_flag=1;
float use_2nd_texture_flag=0;
void loadUniforms(GLuint shader_programme)
{	
/*---------------------------SET RENDERING DEFAULTS---------------------------*/

	// Choose vertex and fragment shaders to use as well as view and proj matrices.
	int model_mat_location = glGetUniformLocation (shader_programme, "model_mat");
	int view_mat_location  = glGetUniformLocation (shader_programme, "view_mat");
	int proj_mat_location  = glGetUniformLocation (shader_programme, "proj_mat");
	
	glUniformMatrix4fv (view_mat_location, 1, GL_FALSE, view_mat.m);
	glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, proj_mat.m);
	glUniformMatrix4fv (model_mat_location, 1, GL_FALSE, model_mat.m);
	
	// WRITE CODE TO LOAD OTHER UNIFORM VARIABLES LIKE FLAGS FOR ENABLING OR DISABLING CERTAIN FUNCTIONALITIES

	int light_position_location = glGetUniformLocation (shader_programme, "light_position");
	int shininess_location = glGetUniformLocation (shader_programme, "shininess");
	int diffuse_effect_flag_location = glGetUniformLocation (shader_programme, "diffuse_effect_flag");
	int specular_effect_flag_location = glGetUniformLocation (shader_programme, "specular_effect_flag");
	int use_texture_flag_location = glGetUniformLocation (shader_programme, "use_texture_flag");
	int use_2nd_texture_flag_location = glGetUniformLocation (shader_programme, "use_2nd_texture_flag");

	glUniform3f(light_position_location, light_x, light_y, light_z);
	glUniform1f(shininess_location, shininess_value);
	glUniform1f(diffuse_effect_flag_location, diffuse_effect_flag);
	glUniform1f(specular_effect_flag_location, specular_effect_flag);
	glUniform1f(use_texture_flag_location, use_texture_flag);
	glUniform1f(use_2nd_texture_flag_location, use_2nd_texture_flag);
}

void drawSurfaceOfRevolution()
{
	// MODIFY THIS LINE OF CODE APPRORIATELY FOR YOUR SURFACE OF REVOLUTION
	glDrawArrays(GL_TRIANGLES, 0, pointCount);
}
	
void keyboardFunction(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// MODIFY THIS FUNCTION FOR KEYBOARD INTERACTIVITY
	//GLFW Reference Links:
	// Callback Example: https://www.glfw.org/docs/3.3/input_guide.html#input_key
	// List of Keys: https://www.glfw.org/docs/3.3/group__keys.html
	
    if (key == GLFW_KEY_A && action == GLFW_PRESS)
    {
	printf("\nKey 'A' pressed, move light source left.\n");
	light_x -= light_change;
    }
    if (key == GLFW_KEY_S && action == GLFW_PRESS)
    {
	printf("\nKey 'S' pressed, move light source down.\n");
        light_y -= light_change;
    }
    if (key == GLFW_KEY_D && action == GLFW_PRESS)
    {
	printf("\nKey 'D' pressed, move light source right.\n");
	light_x += light_change;
    }
    if (key == GLFW_KEY_W && action == GLFW_PRESS)
    {
	printf("\nKey 'W' pressed, move light source up.\n");
        light_y += light_change;
    }
    if (key == GLFW_KEY_Q && action == GLFW_PRESS)
    {
	printf("\nKey 'Q' pressed, move light source closer to origin.\n");
	light_z -= light_change;
    }
    if (key == GLFW_KEY_E && action == GLFW_PRESS)
    {
	printf("\nKey 'E' pressed, move light source farther to origin.\n");
	light_z += light_change;
    }
    if (key == GLFW_KEY_U && action == GLFW_PRESS)
    {
	printf("\nKey 'U' pressed, decrease specular exponent. \n");
	shininess_value -= shininess_change;
    }
    if (key == GLFW_KEY_I && action == GLFW_PRESS)
    {
	printf("\nKey 'I' pressed, increase specular exponent.\n");
        shininess_value += shininess_change;
    }
    if (key == GLFW_KEY_J && action == GLFW_PRESS)
    {
	printf("\nKey 'J' pressed, enable/disable diffuse shading effect.\n");
	if (diffuse_effect_flag)
	    diffuse_effect_flag = 0;
	else
	    diffuse_effect_flag = 1;
    }
    if (key == GLFW_KEY_K && action == GLFW_PRESS)
    {
	printf("\nKey 'K' pressed, enable/disable specular lighting effect.\n");
        if (specular_effect_flag)
	    specular_effect_flag = 0;
	else
	    specular_effect_flag = 1;
    }
    if (key == GLFW_KEY_L && action == GLFW_PRESS)
    {
	printf("\nKey 'L' pressed, enable/disable 1st texture effect.\n");
        if (use_texture_flag)
	    use_texture_flag = 0;
	else
	    use_texture_flag = 1;
    }
    if (key == GLFW_KEY_H && action == GLFW_PRESS)
    {
	printf("\nKey 'H' pressed, enable/disable 2nd texture as normal effect.\n");
        if (use_2nd_texture_flag)
	    use_2nd_texture_flag = 0;
	else
	    use_2nd_texture_flag = 1;
    }	
    if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
	printf("\nKey 'R' pressed, reset everything.\n");
        light_x=0, light_y=0, light_z=5;
	shininess_value=2;
	diffuse_effect_flag=1;
	specular_effect_flag=1;
	use_texture_flag=1;
	use_2nd_texture_flag=0;
    }
	
    if (GLFW_PRESS == glfwGetKey (g_window, GLFW_KEY_ESCAPE)) {
	// Close window when esacape is pressed		
	glfwSetWindowShouldClose (g_window, 1);
    }

}

void binomialCoeffs(int n, int* C){
    int k, j;
    for (k=0; k<= n; k++){
        // compute n!/(k!(n-k)!)
         C[k] = 1;
         for (j=n; j>=k+1; j--)
             C[k] *= j;
         for (j=n-k; j>=2; j--)
              C[k] /= j;
    }
}

GLfloat* getRevolutionVertices(float* cpX, float* cpY, int nCtrlPts,  int ysteps, int thetasteps){
    float *pairY, *pairR;   // save (y,r)pair
    pairY = new float[ysteps+1]; pairR = new float[ysteps+1];
    // use 5th-order Bezier to calculate (y, r) pair for ysteps times
    int n, nBezCurvePts, *C, k, j;
    float tmpX=0.0, tmpY=0.0;
    float u, bezBlendFcn;
    n = nCtrlPts - 1;
    nBezCurvePts = ysteps;
    // binomial coefficients
    C = new int[nCtrlPts];
    binomialCoeffs(nCtrlPts-1, C);
    for (k=0; k<=nBezCurvePts; k++){
        u = float(k) / nBezCurvePts;
        tmpX=0.0, tmpY=0.0;
        for (j=0; j<nCtrlPts; j++){
            bezBlendFcn = C[j]*pow(u,j)*pow(1-u, n-j);
            tmpX += cpX[j] * bezBlendFcn;
            tmpY += cpY[j] * bezBlendFcn;
        }
        //drawPixel(std::round(tmpX), std::round(tmpY));
        pairY[k] = tmpY; pairR[k] = tmpX;
    }
    delete[] C;
    // after get (y, r) pair, start drawing revolution with 2 for loops
    float theta1, theta2, posx1, posz1, posx2, posz2;
    float delta = 2 * MY_PI / thetasteps;
    GLfloat * myVp = new GLfloat[6*3*ysteps*thetasteps];
    int pntCnt = 0; 
    for (j = 0; j < thetasteps; j++){
        theta1 = j * delta;
        theta2 = theta1 + delta;
        //glColor3f (0.0, 0.5*theta1/MY_PI, 0.0);                 
        for (k = 0; k < ysteps; k++){
            posx1 = pairR[k]*cos(theta1); posz1 = pairR[k]*sin(theta1);
            posx2 = pairR[k]*cos(theta2); posz2 = pairR[k]*sin(theta2);
	    //pntCnt += 18 * j * k;
            //face 1, vertex 1
	    myVp[pntCnt++] = posx1; //x
	    myVp[pntCnt++] = pairY[k+1]; //y
	    myVp[pntCnt++] = posz1; //z
	    //face 1, vertex 2
	    myVp[pntCnt++] = posx2; //x
	    myVp[pntCnt++] = pairY[k+1]; //y
	    myVp[pntCnt++] = posz2; //z
	    //face 1, vertex 3
	    myVp[pntCnt++] = posx1; //x
	    myVp[pntCnt++] = pairY[k]; //y
	    myVp[pntCnt++] = posz1; //z
		
	    //face 2, vertex 1
	    myVp[pntCnt++] = posx1; //x
	    myVp[pntCnt++] = pairY[k]; //y
	    myVp[pntCnt++] = posz1; //z
	    //face 2, vertex 2
	    myVp[pntCnt++] = posx2; //x
	    myVp[pntCnt++] = pairY[k+1]; //y
	    myVp[pntCnt++] = posz2; //z
	    //face 2, vertex 3
	    myVp[pntCnt++] = posx2; //x
	    myVp[pntCnt++] = pairY[k]; //y
	    myVp[pntCnt++] = posz2; //z
        }
    }
    delete[] pairY;
    delete[] pairR;
    return myVp;
}

GLfloat* getRevolutionNormals(float* cpX, float* cpY, int nCtrlPts,  int ysteps, int thetasteps){
    float *pairY, *pairR;   // save (y,r)pair
    pairY = new float[ysteps+1]; pairR = new float[ysteps+1];
    // use 5th-order Bezier to calculate (y, r) pair for ysteps times
    int n, nBezCurvePts, *C, k, j;
    float tmpX=0.0, tmpY=0.0;
    float u, bezBlendFcn;
    n = nCtrlPts - 1;
    nBezCurvePts = ysteps;
    // binomial coefficients
    C = new int[nCtrlPts];
    binomialCoeffs(nCtrlPts-1, C);
    for (k=0; k<=nBezCurvePts; k++){
        u = float(k) / nBezCurvePts;
        tmpX=0.0, tmpY=0.0;
        for (j=0; j<nCtrlPts; j++){
            bezBlendFcn = C[j]*pow(u,j)*pow(1-u, n-j);
            tmpX += cpX[j] * bezBlendFcn;
            tmpY += cpY[j] * bezBlendFcn;
        }
        //drawPixel(std::round(tmpX), std::round(tmpY));
        pairY[k] = tmpY; pairR[k] = tmpX;
    }
    delete[] C;
    // after get (y, r) pair, start drawing revolution with 2 for loops
    float theta1, posx1, posz1;
    float delta = 2 * MY_PI / thetasteps;
    int dim1 = thetasteps + 1, dim2 = ysteps + 1, dim3 = 3;
    //float *allPointsIndices = new float[thetasteps+1][ysteps+1][3];
    float *allPointsIndices = new float[dim1*dim2*dim3];
    for (j = 0; j <= thetasteps; j++){
        theta1 = j * delta;
        //glColor3f (0.0, 0.5*theta1/MY_PI, 0.0);                 
        for (k = 0; k <= ysteps; k++){
            posx1 = pairR[k]*cos(theta1); posz1 = pairR[k]*sin(theta1);
            allPointsIndices[j*dim2*dim3+k*dim3+0] = posx1;
	    allPointsIndices[j*dim2*dim3+k*dim3+1] = pairY[k];
            allPointsIndices[j*dim2*dim3+k*dim3+2] = posz1;
	}
    }
    delete[] pairY;
    delete[] pairR;
    float *allPointsNormals = new float[dim1*dim2*dim3];
    vec3 a(0,0,0), b(0,0,0), c(0,0,0), d(0,0,0), normal(0,0,0), zeroVec3(0,0,0); 
    for (j = 0; j <= thetasteps; j++){
        for (k = 0; k <= ysteps; k++){
	    // vector a point to left, b point to bottom, c point to right, d point to top
	    if (j - 1 < 0){
	        a.v[0] = allPointsIndices[(thetasteps-1)*dim2*dim3+k*dim3+0] - allPointsIndices[j*dim2*dim3+k*dim3+0];
	        a.v[1] = allPointsIndices[(thetasteps-1)*dim2*dim3+k*dim3+1] - allPointsIndices[j*dim2*dim3+k*dim3+1];
	        a.v[2] = allPointsIndices[(thetasteps-1)*dim2*dim3+k*dim3+2] - allPointsIndices[j*dim2*dim3+k*dim3+2];
	    }
	    else{
	        a.v[0] = allPointsIndices[(j-1)*dim2*dim3+k*dim3+0] - allPointsIndices[j*dim2*dim3+k*dim3+0];
	        a.v[1] = allPointsIndices[(j-1)*dim2*dim3+k*dim3+1] - allPointsIndices[j*dim2*dim3+k*dim3+1];
	        a.v[2] = allPointsIndices[(j-1)*dim2*dim3+k*dim3+2] - allPointsIndices[j*dim2*dim3+k*dim3+2];
	    }
	    if (j + 1 > thetasteps){
		c.v[0] = allPointsIndices[1*dim2*dim3+k*dim3+0] - allPointsIndices[j*dim2*dim3+k*dim3+0];
		c.v[1] = allPointsIndices[1*dim2*dim3+k*dim3+1] - allPointsIndices[j*dim2*dim3+k*dim3+1];
		c.v[2] = allPointsIndices[1*dim2*dim3+k*dim3+2] - allPointsIndices[j*dim2*dim3+k*dim3+2];
	    }
	    else{
	    	c.v[0] = allPointsIndices[(j+1)*dim2*dim3+k*dim3+0] - allPointsIndices[j*dim2*dim3+k*dim3+0];
	    	c.v[1] = allPointsIndices[(j+1)*dim2*dim3+k*dim3+1] - allPointsIndices[j*dim2*dim3+k*dim3+1];
	    	c.v[2] = allPointsIndices[(j+1)*dim2*dim3+k*dim3+2] - allPointsIndices[j*dim2*dim3+k*dim3+2];
	    }
	    if (k - 1 < 0)
		d = zeroVec3;
	    else{
		d.v[0] = allPointsIndices[j*dim2*dim3+(k-1)*dim3+0] - allPointsIndices[j*dim2*dim3+k*dim3+0];
		d.v[1] = allPointsIndices[j*dim2*dim3+(k-1)*dim3+1] - allPointsIndices[j*dim2*dim3+k*dim3+1];
		d.v[2] = allPointsIndices[j*dim2*dim3+(k-1)*dim3+2] - allPointsIndices[j*dim2*dim3+k*dim3+2];
	    }
	    if (k + 1 > ysteps)
		b = zeroVec3;
	    else{
		b.v[0] = allPointsIndices[j*dim2*dim3+(k+1)*dim3+0] - allPointsIndices[j*dim2*dim3+k*dim3+0];
		b.v[1] = allPointsIndices[j*dim2*dim3+(k+1)*dim3+1] - allPointsIndices[j*dim2*dim3+k*dim3+1];
		b.v[2] = allPointsIndices[j*dim2*dim3+(k+1)*dim3+2] - allPointsIndices[j*dim2*dim3+k*dim3+2];
	    }
	    normal = normalise(normalise(cross(a, b))+normalise(cross(b, c))+normalise(cross(c, d))+normalise(cross(d, a)));
	    allPointsNormals[j*dim2*dim3+k*dim3+0] = normal.v[0];            
	    allPointsNormals[j*dim2*dim3+k*dim3+1] = normal.v[1];            
	    allPointsNormals[j*dim2*dim3+k*dim3+2] = normal.v[2];
	}
    }
    delete [] allPointsIndices;
    GLfloat * myVp = new GLfloat[6*3*ysteps*thetasteps];
    int pntCnt = 0;
    for (j = 0; j < thetasteps; j++){
        for (k = 0; k < ysteps; k++){
	    //pntCnt += 18 * j * k;
            //face 1, vertex 1
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+(k+1)*dim3+0]; //x
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+(k+1)*dim3+1]; //y
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+(k+1)*dim3+2]; //z
	    //face 1, vertex 2
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+(k+1)*dim3+0]; //x
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+(k+1)*dim3+1]; //y
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+(k+1)*dim3+2]; //z
	    //face 1, vertex 3
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+k*dim3+0]; //x
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+k*dim3+1]; //y
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+k*dim3+2]; //z
		
	    //face 2, vertex 1
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+k*dim3+0]; //x
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+k*dim3+1]; //y
	    myVp[pntCnt++] = allPointsNormals[j*dim2*dim3+k*dim3+2]; //z
	    //face 2, vertex 2
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+(k+1)*dim3+0]; //x
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+(k+1)*dim3+1]; //y
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+(k+1)*dim3+2]; //z
	    //face 2, vertex 3
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+k*dim3+0]; //x
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+k*dim3+1]; //y
	    myVp[pntCnt++] = allPointsNormals[(j+1)*dim2*dim3+k*dim3+2]; //z
        }
    }
    delete [] allPointsNormals;
    return myVp;
}

GLfloat* getRevolutionTexCoords(int ysteps, int thetasteps){
    GLfloat * myVp = new GLfloat[6*2*ysteps*thetasteps];
    int pntCnt = 0, k, j; 
    for (j = 0; j < thetasteps; j++){
        //glColor3f (0.0, 0.5*theta1/MY_PI, 0.0);                 
        for (k = 0; k < ysteps; k++){

            //face 1, vertex 1
	    myVp[pntCnt++] = 1.0 * j / thetasteps; //x
	    myVp[pntCnt++] = 1.0 * (k+1) / ysteps; //y
	    
	    //face 1, vertex 2
	    myVp[pntCnt++] = 1.0 * (j+1) / thetasteps; //x
	    myVp[pntCnt++] = 1.0 * (k+1) / ysteps; //y
	   
	    //face 1, vertex 3
	    myVp[pntCnt++] = 1.0 * j / thetasteps; //x
	    myVp[pntCnt++] = 1.0 * k / ysteps; //y
	   	
	    //face 2, vertex 1
	    myVp[pntCnt++] = 1.0 * j / thetasteps; //x
	    myVp[pntCnt++] = 1.0 * k / ysteps; //y
	
	    //face 2, vertex 2
	    myVp[pntCnt++] = 1.0 * (j+1) / thetasteps; //x
	    myVp[pntCnt++] = 1.0 * (k+1) / ysteps; //y

	    //face 2, vertex 3
	    myVp[pntCnt++] = 1.0 * (j+1) / thetasteps; //x
	    myVp[pntCnt++] = 1.0 * k / ysteps; //y
        }
    }
    return myVp;
}


