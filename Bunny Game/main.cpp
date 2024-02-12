#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//#include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[3];
int gWidth, gHeight;

GLint modelingMatrixLoc[3];
GLint viewingMatrixLoc[3];
GLint projectionMatrixLoc[3];
GLint eyePosLoc[3];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix;
glm::vec3 eyePos(0, 0, 0);

int activeProgramIndex = 0; // shader set to use



// Bunny
// Set the initial position of the bunny.
glm::vec3 bunnyPosition = glm::vec3(0.0f, 0.0f, 0.0f); // Assuming the bunny starts at the origin.
// Set the initial scene
// Set the camera to follow the bunny's position with some offsets.
float cameraYOffset = 4.5f; // The camera's height above the bunny.
float cameraZOffset = -12.5f; // The camera's distance behind the bunny. Negative value because we're assuming the bunny's forward direction is along the negative Z-axis.
glm::vec3 cameraPos =  glm::vec3(0.0f, 0.0f, 7.0f) + glm::vec3(0.0f, cameraYOffset, cameraZOffset);
glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 200.0f)+  glm::vec3(0.0f, 0.0f, -1.0f); // Camera looks in the direction the bunny is facing.
glm::vec3 upVector = glm::vec3(0.0f, 1.0f, 0.0f); // Up vector.
float speed = 0.4f;
float bunnyHorizontalSpeed = 1.0f;

float ground_start_pos = -10.5f;
float ground_dist = ground_start_pos;
// Set up the modeling matrix for the ground plane
glm::vec3 ground_pos = glm::vec3(0.f, 0.f, ground_dist);
glm::vec3 ground_translate = glm::vec3(0.f, 5.f, ground_dist);
glm::vec3 ground_scale = glm::vec3(9.0f, -200.0f, 0.f);

float speedScaleFactor = 0.4f;
bool canmove = true;

glm::vec3 checkpoint_pos = glm::vec3(0.f, 5.f, 50.f);
glm::vec3 checkpoint_translate = glm::vec3(0.f, 5.f, 50.f);
glm::vec3 checkpoint_scale = glm::vec3(1.0f, 1.0f, 4.0f);

bool set_checkpoint = true; // initially set checkpoint


glm::vec3 bunny_rotation = glm::vec3(0.f, 1.f, 0.f);

// Bunny hitbox
float bunny_hitbox_x = 0.75f;
float bunny_hitbox_y = 0.75f;
float bunny_hitbox_z = 0.75f;


// box hitbox
float box_hitbox_x = 1.0f;
float box_hitbox_y = 1.0f;
float box_hitbox_z = 1.0f;


bool bunny_rotates = false;
bool bunny_faints = false;


float bunny_angle = 90.0f;
float bunny_angle_z = 0.0f;

typedef struct Cube {
	glm::mat4 modelMat;
	int programIndex;
	bool is_yellow;
	glm::vec3 pos;
} Cube;

Cube cubes[3];

struct Vertex
{
	Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};

struct Texture
{
	Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
	GLfloat u, v;
};

struct Normal
{
	Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
	GLuint vIndex[3], tIndex[3], nIndex[3];
};

/*vector<Vertex> gVertices;
vector<Texture> gTextures;
vector<Normal> gNormals;
vector<Face> gFaces;

GLuint gVertexAttribBuffer, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes; */

struct ObjectData
{
	vector<Vertex> vertices;
	vector<Texture> textures;
	vector<Normal> normals;
	vector<Face> faces;
	GLuint vao;
	GLuint vertexAttribBuffer, indexBuffer;
	GLint inVertexLoc, inNormalLoc;
	int vertexDataSizeInBytes, normalDataSizeInBytes;
};

ObjectData bunnyData, quadData, cubeData;

bool ParseObj(const string& fileName, ObjectData& objData)
{
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			stringstream str(curLine);
			GLfloat c1, c2, c3;
			GLuint index[9];
			string tmp;

			if (curLine.length() >= 2)
			{
				if (curLine[0] == 'v')
				{
					if (curLine[1] == 't') // texture
					{
						str >> tmp; // consume "vt"
						str >> c1 >> c2;
						objData.textures.push_back(Texture(c1, c2));
					}
					else if (curLine[1] == 'n') // normal
					{
						str >> tmp; // consume "vn"
						str >> c1 >> c2 >> c3;
						objData.normals.push_back(Normal(c1, c2, c3));
					}
					else // vertex
					{
						str >> tmp; // consume "v"
						str >> c1 >> c2 >> c3;
						objData.vertices.push_back(Vertex(c1, c2, c3));
					}
				}
				else if (curLine[0] == 'f') // face
				{
					str >> tmp; // consume "f"
					char c;
					int vIndex[3], nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0];
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1];
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2];

					assert(vIndex[0] == nIndex[0] &&
						vIndex[1] == nIndex[1] &&
						vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

					objData.faces.push_back(Face(vIndex, tIndex, nIndex));
				}
				else
				{
					cout << "Ignoring unidentified line in obj file: " << curLine << endl;
				}
			}

			//data += curLine;
			if (!myfile.eof())
			{
				//data += "\n";
			}
		}

		myfile.close();
	}
	else
	{
		return false;
	}
	return true;
}

bool ReadDataFromFile(
	const string& fileName, ///< [in]  Name of the shader file
	string& data)     ///< [out] The contents of the file
{
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			data += curLine;
			if (!myfile.eof())
			{
				data += "\n";
			}
		}

		myfile.close();
	}
	else
	{
		return false;
	}

	return true;
}

GLuint createVS(const char* shaderName)
{
	string shaderSource;

	string filename(shaderName);
	if (!ReadDataFromFile(filename, shaderSource))
	{
		cout << "Cannot find file name: " + filename << endl;
		exit(-1);
	}

	GLint length = shaderSource.length();
	const GLchar* shader = (const GLchar*)shaderSource.c_str();

	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vs, 1, &shader, &length);
	glCompileShader(vs);

	char output[1024] = { 0 };
	glGetShaderInfoLog(vs, 1024, &length, output);
	printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
	string shaderSource;

	string filename(shaderName);
	if (!ReadDataFromFile(filename, shaderSource))
	{
		cout << "Cannot find file name: " + filename << endl;
		exit(-1);
	}

	GLint length = shaderSource.length();
	const GLchar* shader = (const GLchar*)shaderSource.c_str();

	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fs, 1, &shader, &length);
	glCompileShader(fs);

	char output[1024] = { 0 };
	glGetShaderInfoLog(fs, 1024, &length, output);
	printf("FS compile log: %s\n", output);

	return fs;
}

void initShaders()
{
	// Create the programs

	gProgram[0] = glCreateProgram();
	gProgram[1] = glCreateProgram();
	gProgram[2] = glCreateProgram();

	// Create the shaders for both programs

	GLuint vs1 = createVS("bunnyvert.glsl");
	GLuint fs1 = createFS("bunnyfrag.glsl");

	GLuint vs2 = createVS("quadvert.glsl");
	GLuint fs2 = createFS("quadfrag.glsl");

	GLuint vs3 = createVS("cubevert.glsl");
	GLuint fs3 = createFS("cubefrag.glsl");
	// Attach the shaders to the programs

	glAttachShader(gProgram[0], vs1);
	glAttachShader(gProgram[0], fs1);

	glAttachShader(gProgram[1], vs2);
	glAttachShader(gProgram[1], fs2);

	glAttachShader(gProgram[2], vs3);
	glAttachShader(gProgram[2], fs3);
	// Link the programs

	glLinkProgram(gProgram[0]);
	GLint status;
	glGetProgramiv(gProgram[0], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}

	glLinkProgram(gProgram[1]);
	glGetProgramiv(gProgram[1], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}

	glLinkProgram(gProgram[2]);
	glGetProgramiv(gProgram[2], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}
	// Get the locations of the uniform variables from both programs

	for (int i = 0; i < 3; ++i)
	{
		modelingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "modelingMatrix");
		viewingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "viewingMatrix");
		projectionMatrixLoc[i] = glGetUniformLocation(gProgram[i], "projectionMatrix");
		eyePosLoc[i] = glGetUniformLocation(gProgram[i], "eyePos");
	}
}

void initVBO(ObjectData& objData) {
    glGenVertexArrays(1, &objData.vao);
    glBindVertexArray(objData.vao);

    glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &objData.vertexAttribBuffer);
	glGenBuffers(1, &objData.indexBuffer);

	assert(objData.vertexAttribBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, objData.vertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, objData.indexBuffer);

	objData.vertexDataSizeInBytes = objData.vertices.size() * 3 * sizeof(GLfloat);
	objData.normalDataSizeInBytes = objData.normals.size() * 3 * sizeof(GLfloat);
	int indexDataSizeInBytes = objData.faces.size() * 3 * sizeof(GLuint);
	
	GLfloat* vertexData = new GLfloat[objData.vertexDataSizeInBytes];
	GLfloat* normalData = new GLfloat[objData.normalDataSizeInBytes];
	GLuint* indexData = new GLuint[indexDataSizeInBytes];

	float minX = 1e6, maxX = -1e6;
	float minY = 1e6, maxY = -1e6;
	float minZ = 1e6, maxZ = -1e6;

	for (int i = 0; i < objData.vertices.size(); ++i)
	{
		vertexData[3 * i] = objData.vertices[i].x;
		vertexData[3 * i + 1] = objData.vertices[i].y;
		vertexData[3 * i + 2] = objData.vertices[i].z;

		minX = min(minX, objData.vertices[i].x);
		maxX = max(maxX, objData.vertices[i].x);
		minY = min(minY, objData.vertices[i].y);
		maxY = max(maxY, objData.vertices[i].y);
		minZ = min(minZ, objData.vertices[i].z);
		maxZ = max(maxZ, objData.vertices[i].z);
	}

	std::cout << "minX = " << minX << std::endl;
	std::cout << "maxX = " << maxX << std::endl;
	std::cout << "minY = " << minY << std::endl;
	std::cout << "maxY = " << maxY << std::endl;
	std::cout << "minZ = " << minZ << std::endl;
	std::cout << "maxZ = " << maxZ << std::endl;

	for (int i = 0; i < objData.normals.size(); ++i)
	{
		normalData[3 * i] = objData.normals[i].x;
		normalData[3 * i + 1] = objData.normals[i].y;
		normalData[3 * i + 2] = objData.normals[i].z;
	}

	for (int i = 0; i < objData.faces.size(); ++i)
	{
		indexData[3 * i] = objData.faces[i].vIndex[0];
		indexData[3 * i + 1] = objData.faces[i].vIndex[1];
		indexData[3 * i + 2] = objData.faces[i].vIndex[2];
	}


	glBufferData(GL_ARRAY_BUFFER, objData.vertexDataSizeInBytes + objData.normalDataSizeInBytes, NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, objData.vertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, objData.vertexDataSizeInBytes, objData.normalDataSizeInBytes, normalData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

	// done copying to GPU memory; can free now from CPU memory
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(objData.vertexDataSizeInBytes));

	// Unbind the VAO
	glBindVertexArray(0);
}


void init()
{
	ParseObj("bunny.obj", bunnyData);
	ParseObj("quad.obj", quadData);
	ParseObj("cube.obj", cubeData);


	glEnable(GL_DEPTH_TEST);
	initShaders();
	initVBO(bunnyData);
	initVBO(quadData);
	initVBO(cubeData);
}

void drawModel(ObjectData &objData)
{
	glBindVertexArray(objData.vao); // Bind the VAO
	//std::cout << objData.faces.size() << std::endl;
    glDrawElements(GL_TRIANGLES, objData.faces.size() * 3, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0); // Unbind the VAO
}

glm::mat4 ground_model(float angle, int axis, glm::vec3 translation, glm::vec3 scale, glm::vec3 position, glm::vec3 camera) {
	
	glm::mat4 translationMatrix = glm::translate(glm::mat4(10.f), position);
	glm::vec3 rotationVector = glm::vec3(0, 0, 0);
	if (axis == 0) // x
		rotationVector.x = 1;
	else if (axis == 1) // y
		rotationVector.y = 1;
	else
		rotationVector.z = 1;
	glm::mat4 rotationMatrix =
		glm::rotate(glm::mat4(1.f), glm::radians(angle), rotationVector);
	glm::mat4 scaleMatrix = glm::scale(glm::mat4(1.f), scale);
	glm::mat4 modelMat = translationMatrix * rotationMatrix * scaleMatrix;

	return modelMat;
}

void hop_bunny(glm::vec3 &position, float speedDelta) {
	// std::cout << "position: " << position.x << ", " << position.y << ", " << position.z << std::endl;
	position.y =  0.5 * sin(glfwGetTime() * 50 * speedDelta) + 1.0f;
}

void forward_bunny(glm::vec3 &position, float speedDelta) {
	float time = glfwGetTime();
	float logTime = log(1 + time); // Adding 1 to avoid log(0) which is undefined

	glm::vec3 direction = glm::vec3(0.0f, 0.0f, -1.0f);

	float scaledLogTime = speedScaleFactor * logTime;

	// Use scaledLogTime instead of time to calculate the distance
	float distance = speedDelta * scaledLogTime;

	// Calculate the new position using the logarithmic time scale
	position -= direction * distance;

	cameraPos -= direction * distance;
	cameraTarget -= direction * distance;

	ground_pos -= direction * distance;

	if (bunny_rotates) {
		bunny_angle += scaledLogTime * 5.0f;
		if (bunny_angle - 90 >= 360.0f) {
			bunny_rotates = false;
			bunny_angle = 90.0f;
		}
	}
	if (bunny_faints) {
		bunny_angle_z += scaledLogTime * 12.0f;
		if (bunny_angle_z - 90 >= 0.0f) {
			bunny_faints = false;
			bunny_angle_z = 90.0f;
		}
	}

}

glm::mat4 cube_model(float angle, int axis, glm::vec3 translation, glm::vec3 scale, glm::vec3 position, glm::vec3 camera) {
	
	glm::mat4 translationMatrix = glm::translate(glm::mat4(10.f), position);
	glm::vec3 rotationVector = glm::vec3(0, 0, 0);
	if (axis == 0) // x
		rotationVector.x = 1;
	else if (axis == 1) // y
		rotationVector.y = 1;
	else
		rotationVector.z = 1;
	glm::mat4 rotationMatrix =
		glm::rotate(glm::mat4(1.f), glm::radians(angle), rotationVector);
	glm::mat4 scaleMatrix = glm::scale(glm::mat4(1.f), scale);
	glm::mat4 modelMat = translationMatrix * rotationMatrix * scaleMatrix;

	return modelMat;
}


void draw_cubes() {
	// From left to right there are 3 cubes
	for (int i = 0; i < 3; i++) {
		activeProgramIndex = cubes[i].programIndex; // todo..
		glm::mat4 modelingMatrix = cubes[i].modelMat;
		glUseProgram(gProgram[activeProgramIndex]);
		GLint cubeColorLoc = glGetUniformLocation(gProgram[activeProgramIndex], "cubeColor");
		if (cubes[i].is_yellow) {
			glUniform1i(cubeColorLoc, 1);
		} else {
			glUniform1i(cubeColorLoc, 0);
		}
		glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
		glUniformMatrix4fv(viewingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(projectionMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniform3fv(eyePosLoc[activeProgramIndex], 1, glm::value_ptr(cameraPos));
		drawModel(cubeData);

	}
}

bool hit(glm::vec3 bunnypos, glm::vec3 boxpos) {
	// Check if bunny is within the cube
	//std::cout << "bunnypos: " << bunnypos.x << ", " << bunnypos.y << ", " << bunnypos.z << std::endl;
	//std::cout << "boxpos: " << boxpos.x << ", " << boxpos.y << ", " << boxpos.z << std::endl;
	if (bunnypos.x + bunny_hitbox_x > boxpos.x - box_hitbox_x && bunnypos.x - bunny_hitbox_x < boxpos.x + box_hitbox_x) {
		return true;
	}
	return false;
} 


void display()
{
	glClearColor(0, 0, 0, 1);
	glClearDepth(1.0f);
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	if (set_checkpoint) {
		// Random integer between 0 and 2
		glm::vec3 checkpoint_x_offset = glm::vec3(6.0f, 0.0f, 0.0f);
		int random = rand() % 3;
		cubes[random].programIndex = 2; // todo: yellow later, set to 3 later
		set_checkpoint = false;
		for (int i = 0; i < 3; i++) {
			if (i != random) {
				cubes[i].programIndex = 2; // todo: red later
				cubes[i].is_yellow = false;
			}
			else {
				cubes[i].is_yellow = true;
			}
			checkpoint_pos = glm::vec3(0,0,bunnyPosition.z) + glm::vec3(0.f, 0.f, 50.f) + checkpoint_x_offset;
			cubes[i].pos = checkpoint_pos;
			// translate the checkpoint
			//std::cout << "temp: " << checkpoint_pos.x << ", " << checkpoint_pos.y << ", " << checkpoint_pos.z << std::endl;
			
			cubes[i].modelMat = cube_model(90, 0, checkpoint_translate, checkpoint_scale, checkpoint_pos, cameraPos);
			checkpoint_x_offset.x -= 6.0f;
		}
	}
	if (bunnyPosition.z - 2.0f > checkpoint_pos.z) {
		set_checkpoint = true;

		// Collision check
		for (int i = 0; i < 3; i++) {
			if (hit(bunnyPosition, cubes[i].pos)) {
				if (cubes[i].is_yellow) {
					// Check if bunny is within the cube
					//std::cout << "happy" << std::endl;
					bunny_rotates = true;

				} else {
					// Check if bunny is within the cube
				
					//std::cout << "game over" << std::endl;
					bunny_faints = true;
					// Wait 1 second before resetting
					speed = 0.0f;
					canmove = false;
					

				}

			}

		}

	}
	// Create the viewing matrix from the camera position, the target position, and the up vector.
	viewingMatrix = glm::lookAt(cameraPos, cameraTarget, upVector);



	activeProgramIndex = 0;
	// Draw the bunny
	glm::mat4 modelingMatrix = glm::mat4(1.0f); // Start with the identity matrix.
	modelingMatrix = glm::translate(modelingMatrix, bunnyPosition);
	modelingMatrix = glm::rotate(modelingMatrix, glm::radians(bunny_angle_z), glm::vec3(0.f, 0.f, 1.f));
	modelingMatrix = glm::rotate(modelingMatrix, glm::radians(bunny_angle), bunny_rotation);
	// Apply the modeling transformations to the bunny here.
	modelingMatrix = glm::scale(modelingMatrix, glm::vec3(0.75f, 0.75f, 0.75f));
	// If the bunny needs to be raised to sit on the ground plane, apply translation here.
	modelingMatrix = glm::translate(modelingMatrix, glm::vec3(0.0f, -0.34f, 0.0f)); // Adjust the Y value as needed to place the bunny on the ground.
	glUseProgram(gProgram[activeProgramIndex]);
	glm::vec3 bunnyColor = glm::vec3(1.0f, 0.9f, 0.2f); // RGB for yellowish color
	GLint objectColorLoc = glGetUniformLocation(activeProgramIndex, "objectColor");
	glUniform3fv(objectColorLoc, 1, glm::value_ptr(bunnyColor));

	glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
	glUniformMatrix4fv(viewingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
	glUniformMatrix4fv(projectionMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
	glUniform3fv(eyePosLoc[activeProgramIndex], 1, glm::value_ptr(cameraPos));
	drawModel(bunnyData);
	

	hop_bunny(bunnyPosition, speed);

	forward_bunny(bunnyPosition, speed);

	// forward_scene(cameraPos, cameraTarget, upVector, speed);


	// Set up the modeling matrix for the ground plane.
	// Scale the ground plane to be 10x10.



	modelingMatrix = ground_model(90, 0, ground_translate, ground_scale, ground_pos, cameraPos);

	activeProgramIndex = 1;
	// Draw the ground checkpoints
	glUseProgram(gProgram[activeProgramIndex]);
	GLint scaleLoc = glGetUniformLocation(gProgram[activeProgramIndex], "scale");
	glUniform1f(scaleLoc, 0.05f);
	GLint offsetLoc = glGetUniformLocation(gProgram[activeProgramIndex], "offset");
	glUniform1f(offsetLoc, 0.0f);
	glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
	glUniformMatrix4fv(viewingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
	glUniformMatrix4fv(projectionMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
	glUniform3fv(eyePosLoc[activeProgramIndex], 1, glm::value_ptr(cameraPos));
	drawModel(quadData);


	// Set up the modeling matrix for the cube
	draw_cubes();


}

void reshape(GLFWwindow* window, int w, int h)
{
	w = w < 1 ? 1 : w;
	h = h < 1 ? 1 : h;

	gWidth = w;
	gHeight = h;

	glViewport(0, 0, w, h);

	// Use perspective projection
	float fovyRad = (float)(90.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, w / (float)h, 0.1f, 200.0f);

	// Assume default camera position and orientation (camera is at
	// (0, 0, 0) with looking at -z direction and its up vector pointing
	// at +y direction)
	// 
	//viewingMatrix = glm::mat4(1);
	viewingMatrix = glm::lookAt(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0) + glm::vec3(0, 0, -1), glm::vec3(0, 1, 0));

}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS || action == GLFW_REPEAT)
	{
		float time;
		float logTime;
		glm::vec3 direction;
		float scaledLogTime;
		float distance;
		switch (key)
		{
		case GLFW_KEY_A:
			if (!canmove) {
				break;
			}
			// Change bunny pos by speed
			time = glfwGetTime();
			logTime = log(1 + time); // Adding 1 to avoid log(0) which is undefined

			direction = glm::vec3(-1.0f, 0.0f, 0.0f);

			scaledLogTime = speedScaleFactor * logTime;

			// Use scaledLogTime instead of time to calculate the distance
			distance = bunnyHorizontalSpeed * scaledLogTime;

			if (bunnyPosition.x - (direction * distance).x < 6.0f) //- ground_pos.x + ground_scale.x )
				bunnyPosition -= direction * distance;
			else
			{
				bunnyPosition.x = 6.0f; //- ground_pos.x + ground_scale.x;
			}
			break;
		case GLFW_KEY_D:
			if (!canmove) {
				break;
			}
			// Change bunny pos by speed
			time = glfwGetTime();
			logTime = log(1 + time); // Adding 1 to avoid log(0) which is undefined

			direction = glm::vec3(1.0f, 0.0f, 0.0f);

			scaledLogTime = speedScaleFactor * logTime;

			// Use scaledLogTime instead of time to calculate the distance
			distance = bunnyHorizontalSpeed * scaledLogTime;


			if (bunnyPosition.x - (direction * distance).x > -6.0f)//ground_pos.x - ground_scale.x )
				bunnyPosition -= direction * distance;	
			else
			{
				bunnyPosition.x = -6.0f; //ground_pos.x - ground_scale.x;
			}		
			break;
		case GLFW_KEY_R:
			// Wait 1 second before resetting
			glfwWaitEventsTimeout(1.0f);


			// Reset bunny position
			bunnyPosition = glm::vec3(0.0f, 0.0f, 0.0f);

			// Reset camera position
			cameraPos = glm::vec3(0.0f, 0.0f, 7.0f) + glm::vec3(0.0f, cameraYOffset, cameraZOffset);

			// Reset camera target
			cameraTarget = glm::vec3(0.0f, 0.0f, 200.0f)+  glm::vec3(0.0f, 0.0f, -1.0f);

			// Reset ground position
			ground_pos = glm::vec3(0.f, 0.f, ground_dist);

			// Reset speed
			speed = 0.4f;
			bunnyHorizontalSpeed = 1.2f;

			// Reset ground scale
			ground_scale = glm::vec3(8.0f, -200.5f, 0.f);

			// Reset time
			glfwSetTime(0.0f);

			// Reset checkpoint
			set_checkpoint = true;
			

			bunny_faints = false;
			bunny_rotates = false;
			bunny_angle = 90.0f;
			bunny_rotation = glm::vec3(0.f, 1.f, 0.f);
			bunny_angle_z = 0.0f;

			canmove = true;

			break;
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;
		}
	}
}

void mainLoop(GLFWwindow* window)
{
	while (!glfwWindowShouldClose(window))
	{
		display();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
	GLFWwindow* window;
	if (!glfwInit())
	{
		exit(-1);
	}

	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this if on MacOS

	int width = 1000, height = 800;
	window = glfwCreateWindow(width, height, "Simple Example", NULL, NULL);

	if (!window)
	{
		glfwTerminate();
		exit(-1);
	}

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	// Initialize GLEW to setup the OpenGL Function pointers
	if (GLEW_OK != glewInit())
	{
		std::cout << "Failed to initialize GLEW" << std::endl;
		return EXIT_FAILURE;
	}

	char rendererInfo[512] = { 0 };
	strcpy(rendererInfo, (const char*)glGetString(GL_RENDERER)); // Use strcpy_s on Windows, strcpy on Linux
	strcat(rendererInfo, " - "); // Use strcpy_s on Windows, strcpy on Linux
	strcat(rendererInfo, (const char*)glGetString(GL_VERSION)); // Use strcpy_s on Windows, strcpy on Linux
	glfwSetWindowTitle(window, rendererInfo);

	init();

	glfwSetKeyCallback(window, keyboard);
	glfwSetWindowSizeCallback(window, reshape);

	reshape(window, width, height); // need to call this once ourselves
	mainLoop(window); // this does not return unless the window is closed

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
