#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <array>
#include <ctime>
#include <random>
#include <glad/glad.h>
#include <GL/glfw3.h>
#include "Particle.h"
#include "Constraint.h"
#include "Geometry.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Shader.h"
#include "Camera.h"


// Properties
GLuint screenWidth = 1024, screenHeight = 720;
int numParticles = 5;
int numParticlesHeight = 1;
int numParticlesWidth = 5;
int clothHeight = 15;
int clothWidth = 10;
float dt = 0.0075f;
int solver = 1;
float bounce = 1.0;
float damping = 100;
float elasticity = 50;
int stiffnes = 5;
glm::vec3 windDir = glm::vec3(0.0f, 0.0f, 0.0f);

//CalCoreModel xCoreModel;

float stDamp = 30;
float stElast = 500;
float shDamp = 15;
float shElast = 150;
float bDamp = 15;
float bElast = 15;


// Function prototypes
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void Do_Movement();
void moveSphere(Sphere& s);
void check_collisions(std::vector<Plane>& Planes, std::vector<Particle>& Particles, Triangle triangle, int pn, std::vector<std::vector<float>>& disant, std::vector<std::vector<float>>& disact, std::vector<Sphere> Spheres);
void createTriangle(std::vector<GLfloat>& vertices, std::vector<GLushort>& indices, Triangle triangle);
void createSphere(std::vector<GLfloat>& vertices, std::vector<GLushort>& indices, float radius, unsigned int rings, unsigned int sectors);
inline void push_indices(std::vector<GLushort>& indices, int sectors, int r, int s);
void calculateStringForces(Particle& p1, Particle& p2);
bool stringDists(std::vector<Particle>& Particles);
void calcWindForce(Particle& p1, Particle& p2, Particle& p3);

void drawSpheres(std::vector<Sphere> Spheres, std::vector<GLushort> sphereIndices, unsigned int sphereVAO, Shader ourShader);
void drawParticle(Particle Particle, std::vector<GLushort> sphereIndices, unsigned int sphereVAO, Shader ourShader);
void drawTriangle(std::vector<GLushort> triIndices, unsigned int triVAO, Shader ourShader);
void drawScenario(unsigned int VAO, Shader ourShader);

glm::vec3 cross(const glm::vec3 &v, const glm::vec3 &f);
float dot(const glm::vec3 &v, const glm::vec3 &f);

glm::vec3 cross(const glm::vec3 &f, const glm::vec3 &v)
{
	return glm::vec3(f[1] * v[2] - f[2] * v[1], f[2] * v[0] - f[0] * v[2], f[0] * v[1] - f[1] * v[0]);
}

float dot(const glm::vec3 &v, const glm::vec3 &f)
{
	return f[0] * v[0] + f[1] * v[1] + f[2] * v[2];
}

// Camera
int N = 0;
Camera camera(glm::vec3(0.0f, 10.0f, N + 25));
bool keys[1024];
GLfloat lastX = 400, lastY = 300;
bool firstMouse = true;

GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;


void main() {

	//xCoreModel.create("dummy");

	//Rand initialization
	srand(static_cast <unsigned> (time(0)));
	std::random_device rd;
	std::mt19937 e2(rd());
	std::uniform_real_distribution<> dist(-500, 500);
	std::uniform_real_distribution<> dist2(-500, 500);
	std::uniform_real_distribution<> dist3(0, 19);
	// Inicializar GLFW
	glfwInit();
	glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(screenWidth, screenHeight, "LearnOpenGL", nullptr, nullptr); // Windowed
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return;
	}
	glfwMakeContextCurrent(window);

	ImGui_ImplGlfwGL2_Init(window, true);
	// Set the required callback functions
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);

	// Options
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return;
	}

	// Viewport dimensions
	glViewport(0, 0, screenWidth, screenHeight);

	// Setup  OpenGL options
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	Shader ourShader("Resources/Shaders/default.vs", "Resources/Shaders/default.frag");

	// Planes -----------------------------------------------------------
	std::vector<Plane> Planes;

	// Floor plane
	glm::vec3 punt(0.0f);
	glm::vec3 normal(0.0f, 1.0f, 0.0f);
	Planes.push_back(Plane(punt, normal));
	// Left plane
	glm::vec3 puntL(-10.0f, 0.0f, 0.0f);
	glm::vec3 normalL(1.0f, 0.0f, 0.0f);
	Planes.push_back(Plane(puntL, normalL));
	// Right plane
	glm::vec3 puntR(10.0f, 0.0f, 0.0f);
	glm::vec3 normalR(-1.0f, 0.0f, 0.0f);
	Planes.push_back(Plane(puntR, normalR));
	// Back plane
	glm::vec3 puntB(0.0f, 0.0f, -10.0f);
	glm::vec3 normalB(0.0f, 0.0f, 1.0f);
	Planes.push_back(Plane(puntB, normalB));
	// Front plane
	glm::vec3 puntFr(0.0f, 0.0f, 10.0f);
	glm::vec3 normalFr(0.0f, 0.0f, -1.0f);
	Planes.push_back(Plane(puntFr, normalFr));
	//--------------------------------------------------------------------

	// Spheres -----------------------------------------------------------
	std::vector<Sphere> Spheres;
	int numSpheres = 0;

	for (int s = 0; s < numSpheres; ++s) {
		glm::vec3 center1(0.0f, 5.0f, 20 - s * 4);
		Spheres.push_back(Sphere(center1, 2));
	}

	glm::vec3 centerS(20.0f, 4.0f, 10.0f);
	Spheres.push_back(Sphere(centerS, 2));
	//--------------------------------------------------------------------

	// Triangles ---------------------------------------------------------
	glm::vec3 punt1(-2.0f, 7.0f, 20.0f);
	glm::vec3 punt2(2.0f, 10.0f, -20.0f);
	glm::vec3 punt3(2.0f, 7.0f, 20.0f);
	Triangle triangle(punt1, punt2, punt3);
	//--------------------------------------------------------------------

	
#pragma region "object_initialization"
	// set up vertex data (and buffer(s)) and configure vertex attributes
	// ------------------------------------------------------------------
	float vertices[] = {
		-0.5f, -0.5f, -0.5f,  0.0f, 0.0f,
		0.5f, -0.5f, -0.5f,  1.0f, 0.0f,
		0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
		0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
		-0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
		-0.5f, -0.5f, -0.5f,  0.0f, 0.0f,

		-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
		0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
		0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
		0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
		-0.5f,  0.5f,  0.5f,  0.0f, 1.0f,
		-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,

		-0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
		-0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
		-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
		-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
		-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
		-0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

		0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
		0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
		0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
		0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
		0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
		0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

		-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
		0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
		0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
		0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
		-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
		-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,

		-0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
		0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
		0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
		0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
		-0.5f,  0.5f,  0.5f,  0.0f, 0.0f,
		-0.5f,  0.5f, -0.5f,  0.0f, 1.0f
	};

	unsigned int VBO, VAO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	// texture coord attribute
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	// Sphere ---------------------------------------------------------
	std::vector<GLfloat> sphereVertices;
	std::vector<GLushort> sphereIndices;

	createSphere(sphereVertices, sphereIndices, 2, 10, 10);

	unsigned int sphereVBO, sphereVAO, sphereEBO;
	glGenVertexArrays(1, &sphereVAO);
	glGenBuffers(1, &sphereVBO);
	glGenBuffers(1, &sphereEBO);

	glBindVertexArray(sphereVAO);

	glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);
	glBufferData(GL_ARRAY_BUFFER, sphereVertices.size() * sizeof(GLfloat), &sphereVertices[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereEBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sphereIndices.size() * sizeof(GLushort), &sphereIndices[0], GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	//------------------------------------------------------------------

	// Triangle --------------------------------------------------------
	std::vector<GLfloat> triVertices;
	std::vector<GLushort> triIndices;

	createTriangle(triVertices, triIndices, triangle);

	unsigned int triVBO, triVAO, triEBO;
	glGenVertexArrays(1, &triVAO);
	glGenBuffers(1, &triVBO);
	glGenBuffers(1, &triEBO);

	glBindVertexArray(triVAO);

	glBindBuffer(GL_ARRAY_BUFFER, triVBO);
	glBufferData(GL_ARRAY_BUFFER, triVertices.size() * sizeof(GLfloat), &triVertices[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triEBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, triIndices.size() * sizeof(GLushort), &triIndices[0], GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	//------------------------------------------------------------------
	// Line ------------------------------------------------------------
	float vs[] = {
		0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f
	};
	unsigned int lVBO, lVAO;
	glGenVertexArrays(1, &lVAO);
	glGenBuffers(1, &lVBO);

	glBindVertexArray(lVAO);
	glBindBuffer(GL_ARRAY_BUFFER, lVBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(vs), vs, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, lVBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	//------------------------------------------------------------------
	dt = 0.0005f;  //simulation time
	windDir = windDir*dt;
	float tini = 0.0f;
	float time = tini;
	float tfinal = 6.0f; //final time of simulation 
						 // particle inicialization
	
	// simulation loop
	int count = 0;
	int nxtPart = 0;
	std::vector<std::vector<float>> disact = std::vector<std::vector<float>>(Planes.size()), disant = std::vector<std::vector<float>>(Planes.size());

	
	std::vector<Particle> Particles;
	std::vector<Constraint> Constraints;
	for (int pn = 0; pn < 10; pn++) {
		float rx = std::floor(dist3(e2));
		float rz = std::floor(dist3(e2));
		float gx = std::floor(dist3(e2));
		float gz = std::floor(dist3(e2));
		float vx = std::floor(dist(e2)) / 100;
		float vz = std::floor(dist(e2)) / 100;
		Particles.push_back(Particle(rx-9.5f, 0.25f, rz-9.5f)); //initial position of the particle

		Particles.at(pn).setLifetime(10000000.0f);
		std::cout << "Lifetime =" << Particles[pn].getLifetime() << std::endl;
		//Particles.at(pn).setBouncing(0.8f);
		//Particles.at(pn).addForce(0, -9.8f, 0);
		Particles.at(pn).setVelocity(vx, 0.0f, vz);
		glm::vec3 v = Particles.at(pn).getVelocity();
		Particles.at(pn).max_speed =  sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
		Particles.at(pn).setPreviousPosition(Particles.at(pn).getCurrentPosition() - Particles.at(pn).getVelocity()*dt);
		Particles.at(pn).aStar(glm::vec2(int(rx), int(rz)), glm::vec2(int(gx), int(gz)));
		//	p.setFixed(true);

		for (int np = 0; np < Planes.size(); np++) {
			disact.at(np).push_back(Planes.at(np).distPoint2Plane(Particles.at(pn).getCurrentPosition()));
			disant.at(np).push_back(Planes.at(np).distPoint2Plane(Particles.at(pn).getCurrentPosition()));
		};
	}

	/*for (int y = 0; y < numParticlesWidth; y++) {
		for (int x= 0; x < numParticlesHeight; x++) {
			int pn = x + numParticlesHeight*y;
			float rx = std::floor(dist(e2)) / 100;
			float rz = std::floor(dist(e2)) / 100;
			Particles.push_back(Particle(clothWidth * (x / (float)numParticlesWidth) - 5, 1.0f , -clothHeight * (y / (float)numParticlesHeight)-2)); //initial position of the particle

			Particles.at(pn).setLifetime(10000000.0f);
			//std::cout << "Lifetime =" << Particles[pn].getLifetime() << std::endl;
			Particles.at(pn).setBouncing(1.0f);
			//Particles.at(pn).addForce(0, -9.8f, 0);
			Particles.at(pn).setVelocity(rx, 0.0f, rz);
			Particles.at(pn).setPreviousPosition(Particles.at(pn).getCurrentPosition() - Particles.at(pn).getVelocity()*dt);
			//	p.setFixed(true);

			for (int np = 0; np < Planes.size(); np++) {
				disact.at(np).push_back(Planes.at(np).distPoint2Plane(Particles.at(pn).getCurrentPosition()));
				disant.at(np).push_back(Planes.at(np).distPoint2Plane(Particles.at(pn).getCurrentPosition()));
			};
		}
	}*/



	//Particles.at(0).setFixed(true);

	bool show_test_window = true;
	bool show_another_window = false;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	while (!glfwWindowShouldClose(window) || time <= tfinal)
	{

		nxtPart++;

		// Set frame time
		GLfloat currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		// Check and call events
		glfwPollEvents();
		Do_Movement();
		moveSphere(Spheres.at(Spheres.size()-1));

		// render
		// ------
		glClearColor(0.4f, 0.4f, 0.4f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // also clear the depth buffer now!

		// activate shader
		ourShader.Use();

		// pass projection matrix to shader (note that in this case it could change every frame)
		// camera/view transformation
		glm::mat4 view = camera.GetViewMatrix();
		glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)screenWidth / (float)screenHeight, 0.1f, 100.0f);
		glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "view"), 1, GL_FALSE, glm::value_ptr(view));
		glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
		glUniform3f(glGetUniformLocation(ourShader.Program, "cameraPos"), camera.Position.x, camera.Position.y, camera.Position.z);

		// Draw Walls
		drawScenario(VAO, ourShader);
		// Draw Spheres
		//drawSpheres(Spheres, sphereIndices, sphereVAO, ourShader);
		// Draw triangle
		//drawTriangle(triIndices, triVAO, ourShader);

		// Draw particlesW
		for (int pn = 0; pn < Particles.size(); pn++) {
			if (Particles.at(pn).getLifetime() > 0) {

				drawParticle(Particles.at(pn), sphereIndices, sphereVAO, ourShader);

				if (Particles.at(pn).road.empty() && Particles.at(pn).getDistanceToWaypoint() <= 0.15) {

					float gx = std::floor(dist3(e2));
					float gz = std::floor(dist3(e2));

					if (gx == 10) gx = 9;

					glm::vec3 st = Particles.at(pn).getWaypoint();
					Particles.at(pn).aStar(glm::vec2(int(st[0] + 9.5), int(st[2] + 9.5)), glm::vec2(int(gx), int(gz)));

					std::cout << Particles.at(pn).road[0][0] << std::endl;

				}
				else if (Particles.at(pn).getDistanceToWaypoint() <= 0.025)  Particles.at(pn).setNextWaypoint();


				if (solver == 0) Particles.at(pn).updateParticle(dt, Particle::UpdateMethod::EulerOrig);
				else if (solver == 1) Particles.at(pn).updateParticle(dt, Particle::UpdateMethod::EulerSemi);
				else if (solver == 2) Particles.at(pn).updateParticle(dt, Particle::UpdateMethod::Verlet);		

				if (dt > 0) Particles.at(pn).setLifetime(Particles.at(pn).getLifetime() - dt);
				else Particles.at(pn).setLifetime(Particles.at(pn).getLifetime() + dt);

				check_collisions(Planes, Particles, triangle, pn, disant, disact, Spheres);

				for (int pn2 = 0; pn2 < Particles.size(); pn2++) {
					if (pn != pn2) {
						glm::vec3 fPos1 = Particles.at(pn).getCurrentPosition() + Particles.at(pn).getVelocity() * 200.0f * dt;
						glm::vec3 fPos2 = Particles.at(pn2).getCurrentPosition() + Particles.at(pn2).getVelocity() * 200.0f * dt;

						//glm::vec3 distVec = fPos2 - fPos1;
						//float dist = distVec.length();

						glm::vec3 p1_to_p2 = fPos2 - fPos1;
						float norm = sqrt(p1_to_p2.x * p1_to_p2.x + p1_to_p2.y * p1_to_p2.y + p1_to_p2.z * p1_to_p2.z); // current distance between p1 and p2
						p1_to_p2 = p1_to_p2 / norm;

						if (norm < 0.2f) {
							Particles.at(pn).avoidance = glm::vec3(p1_to_p2.z, p1_to_p2.y, -p1_to_p2.x);
						}
					}
				}


			}
			else {

				Particles.erase(Particles.begin() + pn);
				for (int np = 0; np < Planes.size(); np++) {
					disant.at(np).erase(disant.at(np).begin() + pn*Planes.size());
					disact.at(np).erase(disact.at(np).begin() + pn*Planes.size());
				}
				
				pn = pn - 1;
			}
		}
		


		if(dt > 0) time = time + dt; //increase time counter
		else time = time - dt;

		glfwSwapBuffers(window);
	}
	//system("PAUSE");
	ImGui_ImplGlfwGL2_Shutdown();
	glfwTerminate();
	return;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void check_collisions(std::vector<Plane>& Planes, std::vector<Particle>& Particles, Triangle triangle, int pn, std::vector<std::vector<float>>& disant, std::vector<std::vector<float>>& disact, std::vector<Sphere> Spheres)
{
	for (int np = 0; np < Planes.size(); np++) {
		disant.at(np).at(pn) = disact.at(np).at(pn);
		disact.at(np).at(pn) = Planes.at(np).distPoint2Plane(Particles.at(pn).getCurrentPosition());

		if (disant.at(np).at(pn)*disact.at(np).at(pn) < 0.0f) {
			/*glm::vec3 v = reflect(Particles.at(pn).getVelocity(), Planes.at(np).normal);
			float l = (v + Particles.at(pn).getVelocity()).length() * 0.5f;
			Particles.at(pn).setVelocity((v - (l * normalize(Planes.at(np).normal) * (1 - bounce))) * 0.8f);*/
			//Particles.at(pn).setPosition(Particles.at(pn).getPreviousPosition() + Particles.at(pn).getVelocity()*dt);
			//Particles.at(pn).setPosition(Particles.at(pn).getCurrentPosition() - (1.0f + Particles.at(pn).getBouncing()) * (Planes.at(np).normal * Particles.at(pn).getCurrentPosition() + disact.at(np).at(pn)) * Planes.at(np).normal);
			//Particles.at(pn).setVelocity(Particles.at(pn).getVelocity() - (1.0f + Particles.at(pn).getBouncing()) * (Planes.at(np).normal * Particles.at(pn).getVelocity()) * Planes.at(np).normal);
			Particles.at(pn).setVelocity(bounce * reflect(Particles.at(pn).getVelocity(), Planes.at(np).normal));
			glm::vec3 midpoint = (Particles.at(pn).getPreviousPosition() + Particles.at(pn).getCurrentPosition()) / 2.0f;
			Particles.at(pn).setPreviousPosition(midpoint);
			Particles.at(pn).setPosition(midpoint + Particles.at(pn).getVelocity()*dt);
			disact.at(np).at(pn) = disant.at(np).at(pn);
		}
	}

	for (int s = 0; s < Spheres.size(); ++s) {
		if (Spheres.at(s).isInside(Particles.at(pn).getCurrentPosition())) {
			Particles.at(pn).setVelocity(bounce * reflect(Particles.at(pn).getVelocity(), normalize(Particles.at(pn).getCurrentPosition() - Spheres.at(s).center)));
			//Particles.at(pn).setPosition(Particles.at(pn).getPreviousPosition());
			glm::vec3 midpoint = (Particles.at(pn).getPreviousPosition() + Particles.at(pn).getCurrentPosition()) / 2.0f;
			Particles.at(pn).setPreviousPosition(midpoint);
			Particles.at(pn).setPosition(midpoint + Particles.at(pn).getVelocity()*dt);
		}
	}
	glm::vec3 triIntersection;
	if (triangle.intersecSegment(Particles.at(pn).getCurrentPosition(), Particles.at(pn).getPreviousPosition(), triIntersection)) {

		Particles.at(pn).setVelocity(bounce * reflect(Particles.at(pn).getVelocity(), triangle.normal));
		glm::vec3 midpoint = (Particles.at(pn).getPreviousPosition() + Particles.at(pn).getCurrentPosition()) / 2.0f;
		Particles.at(pn).setPreviousPosition(midpoint);
		Particles.at(pn).setPosition(midpoint + Particles.at(pn).getVelocity()*dt);
	}
}

void calculateStringForces(Particle& par1, Particle& par2) {

	glm::vec3 p1 = par1.getCurrentPosition();
	glm::vec3 p2 = par2.getCurrentPosition();
	glm::vec3 v1 = par1.getVelocity();
	glm::vec3 v2 = par2.getVelocity();
	float l = (p2 - p1).length();
	glm::vec3 n = (p2 - p1) / l;

	float fe = elasticity * (l - 1.00f);
	glm::vec3 fd = damping * (v2 - v1) * n;
	glm::vec3 f = (fe + fd) * n;

	par1.addStringfForce(f);
	par2.addStringfForce(-f);
}

void calcWindForce(Particle& p1, Particle& p2, Particle& p3)
{
	glm::vec3 pos1 = p1.getCurrentPosition();
	glm::vec3 pos2 = p2.getCurrentPosition();
	glm::vec3 pos3 = p3.getCurrentPosition();

	glm::vec3 v1 = pos2 - pos1;
	glm::vec3 v2 = pos3 - pos1;

	glm::vec3 normal = cross(v1, v2);

	glm::vec3 d = normal/sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
	glm::vec3 force = normal*(dot(d, windDir));

	p1.addStringfForce(force);
	p2.addStringfForce(force);
	p3.addStringfForce(force);
}

bool stringDists(std::vector<Particle>& Particles) {

	std::vector<float> aux;
	for (int pn = numParticles; pn < numParticles + 25 - 1; pn++) {
		if (pn % 5 == 4) {
			aux.push_back((Particles.at(pn).getCurrentPosition() - Particles.at(pn + 5).getCurrentPosition()).length());
			//aux.push_back((Particles.at(pn).getCurrentPosition() - Particles.at(pn + 4).getCurrentPosition()).length());
		}
		else if ((pn - numParticles - 25) >= -5)
			aux.push_back((Particles.at(pn).getCurrentPosition() - Particles.at(pn + 1).getCurrentPosition()).length());
		else {
			aux.push_back((Particles.at(pn).getCurrentPosition() - Particles.at(pn + 1).getCurrentPosition()).length());
			aux.push_back((Particles.at(pn).getCurrentPosition() - Particles.at(pn + 5).getCurrentPosition()).length());
			//aux.push_back((Particles.at(pn).getCurrentPosition() - Particles.at(pn + 6).getCurrentPosition()).length());
			//if (pn % 5 >= 1)
				//aux.push_back((Particles.at(pn).getCurrentPosition() - Particles.at(pn + 4).getCurrentPosition()).length());
		}
	}

	float max = *max_element(aux.begin(), aux.end());
	float min = *min_element(aux.begin(), aux.end());

	if ((max > 3) || (min <= 2))
		return true;
	return false;
}

#pragma region "Draw elements"
void drawParticle(Particle Particle, std::vector<GLushort> sphereIndices, unsigned int sphereVAO, Shader ourShader)
{
	glBindVertexArray(sphereVAO);
	glm::mat4 model;
	model = glm::translate(model, Particle.getCurrentPosition());
	model = glm::scale(model, glm::vec3(0.1f, 0.1f, 0.1f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(model));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 1.0f, 0.0f, 0.0f, 1.0f);

	glDrawElements(GL_TRIANGLES, sphereIndices.size(), GL_UNSIGNED_SHORT, 0);
}

void drawSpheres(std::vector<Sphere> Spheres, std::vector<GLushort> sphereIndices, unsigned int sphereVAO, Shader ourShader)
{
	for (int s = 0; s < Spheres.size(); ++s) {
		glBindVertexArray(sphereVAO);
		glm::mat4 sphereD;
		sphereD = glm::translate(sphereD, Spheres.at(s).center);
		//sphereD = glm::scale(sphereD, glm::vec3(10.0f, 10.0f, 10.0f));
		glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(sphereD));
		glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.5, 1.0f - s / Spheres.size(), 1.0f - 1.0f*s / Spheres.size(), 0.8f);

		glDrawElements(GL_TRIANGLES, sphereIndices.size(), GL_UNSIGNED_SHORT, 0);
	}
}

void drawTriangle(std::vector<GLushort> triIndices, unsigned int triVAO, Shader ourShader)
{
	glBindVertexArray(triVAO);
	glm::mat4 tri;
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(tri));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.5, 0.4, 1.0, 0.8f);

	glDrawElements(GL_TRIANGLES, triIndices.size(), GL_UNSIGNED_SHORT, 0);
}

void drawScenario(unsigned int VAO, Shader ourShader)
{
	// render Walls
	glBindVertexArray(VAO);

	glm::mat4 floor;
	floor = glm::translate(floor, glm::vec3(0.0, -0.5, 0.0));
	float angle = 90.0f;
	floor = glm::rotate(floor, glm::radians(angle), glm::vec3(1.0f, 0.0f, 0.0f));
	floor = glm::scale(floor, glm::vec3(20.0f, 20.0f, 1.0f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(floor));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 1.0f, 1.0f, 1.0f, 1.0f);

	glDrawArrays(GL_TRIANGLES, 0, 6);

	glm::mat4 wallL;
	wallL = glm::translate(wallL, glm::vec3(-9.5, 1.0, 0.0));
	angle = 90.0f;
	wallL = glm::rotate(wallL, glm::radians(angle), glm::vec3(0.0f, 1.0f, 0.0f));
	wallL = glm::scale(wallL, glm::vec3(20.0f, 2.0f, 1.0f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(wallL));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.7f, 0.7f, 0.7f, 0.8f);

	glDrawArrays(GL_TRIANGLES, 0, 6);

	glm::mat4 wallR;
	wallR = glm::translate(wallR, glm::vec3(9.5, 1.0, 0.0));
	angle = -90.0f;
	wallR = glm::rotate(wallR, glm::radians(angle), glm::vec3(0.0f, 1.0f, 0.0f));
	wallR = glm::scale(wallR, glm::vec3(20.0f, 2.0f, 1.0f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(wallR));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.7f, 0.7f, 0.7f, 0.8f);

	glDrawArrays(GL_TRIANGLES, 0, 6);

	glm::mat4 wallB;
	wallB = glm::translate(wallB, glm::vec3(0.0, 1.0, -9.5));
	wallB = glm::scale(wallB, glm::vec3(20.0f, 2.0f, 1.0f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(wallB));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.75f, 0.75f, 0.75f, 0.8f);

	glDrawArrays(GL_TRIANGLES, 0, 6);

	glm::mat4 wallFr;
	wallFr = glm::translate(wallFr, glm::vec3(0.0, 1.0, 10.5));
	wallFr = glm::scale(wallFr, glm::vec3(20.0f, 2.0f, 1.0f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(wallFr));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.75f, 0.75f, 0.75f, 0.8f);

	glDrawArrays(GL_TRIANGLES, 0, 6);

	glm::mat4 wallUpperMid;
	wallUpperMid = glm::translate(wallUpperMid, glm::vec3(0.5, 0.5, -6.0));
	wallUpperMid = glm::scale(wallUpperMid, glm::vec3(1.0f, 1.0f, 8.0f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(wallUpperMid));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.75f, 0.75f, 0.75f, 0.8f);

	glDrawArrays(GL_TRIANGLES, 0, 36);

	glm::mat4 wallLMid;
	wallLMid = glm::translate(wallLMid, glm::vec3(0.5, 0.5, 6.0));
	wallLMid = glm::scale(wallLMid, glm::vec3(1.0f, 1.0f, 8.0f));
	glUniformMatrix4fv(glGetUniformLocation(ourShader.Program, "model"), 1, GL_FALSE, glm::value_ptr(wallLMid));
	glUniform4f(glGetUniformLocation(ourShader.Program, "FragColor"), 0.75f, 0.75f, 0.75f, 0.8f);

	glDrawArrays(GL_TRIANGLES, 0, 36);

}

#pragma region "Geometry initialitzation"
inline void push_indices(std::vector<GLushort>& indices, int sectors, int r, int s) {
	int curRow = r * sectors;
	int nextRow = (r + 1) * sectors;
	int nextS = (s + 1) % sectors;

	indices.push_back(curRow + s);
	indices.push_back(nextRow + s);
	indices.push_back(nextRow + nextS);

	indices.push_back(curRow + s);
	indices.push_back(nextRow + nextS);
	indices.push_back(curRow + nextS);
}

void createSphere(std::vector<GLfloat>& vertices, std::vector<GLushort>& indices, float radius, unsigned int rings, unsigned int sectors)
{
	float const R = 1. / (float)(rings - 1);
	float const S = 1. / (float)(sectors - 1);

	for (int r = 0; r < rings; ++r) {
		for (int s = 0; s < sectors; ++s) {
			float const y = sin(-M_PI_2 + M_PI * r * R);
			float const x = cos(2 * M_PI * s * S) * sin(M_PI * r * R);
			float const z = sin(2 * M_PI * s * S) * sin(M_PI * r * R);

			//texcoords.push_back(vec2(s*S, r*R));
			vertices.push_back(x * radius);
			vertices.push_back(y * radius);
			vertices.push_back(z * radius);
			if (r < rings - 1)
				push_indices(indices, sectors, r, s);
		}	
	}
}

void createTriangle(std::vector<GLfloat>& vertices, std::vector<GLushort>& indices, Triangle triangle) {
	
	vertices.push_back(triangle.vertex1.x);
	vertices.push_back(triangle.vertex1.y);
	vertices.push_back(triangle.vertex1.z);

	vertices.push_back(triangle.vertex2.x);
	vertices.push_back(triangle.vertex2.y);
	vertices.push_back(triangle.vertex2.z);

	vertices.push_back(triangle.vertex3.x);
	vertices.push_back(triangle.vertex3.y);
	vertices.push_back(triangle.vertex3.z);

	indices.push_back(0);
	indices.push_back(1);
	indices.push_back(2);
}

#pragma region "User input"

// Moves/alters the camera positions based on user input
void Do_Movement()
{
	// Camera controls
	if (keys[GLFW_KEY_W])
		camera.ProcessKeyboard(FORWARD, deltaTime * 5);
	if (keys[GLFW_KEY_S])
		camera.ProcessKeyboard(BACKWARD, deltaTime * 5);
	if (keys[GLFW_KEY_A])
		camera.ProcessKeyboard(LEFT, deltaTime * 5);
	if (keys[GLFW_KEY_D])
		camera.ProcessKeyboard(RIGHT, deltaTime * 5);
	if (keys[GLFW_KEY_R])
		dt += 0.00005f;
	if (keys[GLFW_KEY_F])
		dt -= 0.00005f;
	if (keys[GLFW_KEY_Z]) {
		std::cout << "Solver: Euler Original" << std::endl;
		solver = 0;
	}
	if (keys[GLFW_KEY_X]) {
		std::cout << "Solver: Euler Semi-implicit" << std::endl;
		solver = 1;
	}
	if (keys[GLFW_KEY_C]) {
		std::cout << "Solver: Verlet" << std::endl;
		solver = 2;
	}
	if (keys[GLFW_KEY_B]) {
		if (bounce < 1) bounce += 0.05f*deltaTime;
		std::cout << "Bounce: " << bounce << std::endl;
	}
	if (keys[GLFW_KEY_V]) {
		if (bounce > 0) bounce -= 0.05f*deltaTime;
		std::cout << "Bounce: " << bounce << std::endl;
	}

	if (keys[GLFW_KEY_T]) {
		elasticity += 20.5f*deltaTime;
		std::cout << "Elasticity: " << elasticity << std::endl;
	}
	if (keys[GLFW_KEY_G]) {
		if (elasticity > 0) elasticity -= 20.5f*deltaTime;
		std::cout << "Elasticity: " << elasticity << std::endl;
	}
	if (keys[GLFW_KEY_Y]) {
		damping += 20.5f*deltaTime;
		std::cout << "Damping: " << damping << std::endl;
	}
	if (keys[GLFW_KEY_H]) {
		if (damping > 0) damping -= 20.5f*deltaTime;
		std::cout << "Damping: " << damping << std::endl;
	}
}

void moveSphere(Sphere& s)
{
	if (keys[GLFW_KEY_I])
		s.setPosition(s.center + glm::vec3(0, 2, 0)*dt);
	if (keys[GLFW_KEY_J])
		s.setPosition(s.center + glm::vec3(-2, 0, 0)*dt);
	if (keys[GLFW_KEY_K])
		s.setPosition(s.center + glm::vec3(0, -2, 0)*dt);
	if (keys[GLFW_KEY_L])
		s.setPosition(s.center + glm::vec3(2, 0, 0)*dt);
}

// Callbacks for camera movement
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	if (action == GLFW_PRESS)
		keys[key] = true;
	else if (action == GLFW_RELEASE)
		keys[key] = false;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	GLfloat xoffset = xpos - lastX;
	GLfloat yoffset = lastY - ypos;

	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	camera.ProcessMouseScroll(yoffset);
}
