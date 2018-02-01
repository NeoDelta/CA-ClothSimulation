#include <iostream>
#include <GL/freeglut.h>
#include "Particle.h"
#include "Geometry.h"

// Function Prototypes
void display();

int notmain(int argc, char* argv[]) {

	//  Initialize GLUT and process user parameters
	glutInit(&argc, argv);
	//  Request double buffered true color window with Z-buffer
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	// Create window
	glutCreateWindow("Falling particles");
	//  Enable Z-buffer depth test
	glEnable(GL_DEPTH_TEST);
	// Callback functions
	glutDisplayFunc(display);

	//  Pass control to GLUT for events
	glutMainLoop();

	//  Return to OS
	return 0;
}

void display() {

	//  Clear screen and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();


	//Multi-colored side - FRONT
	glBegin(GL_POLYGON);

	glColor3f(1.0, 0.0, 0.0);     glVertex3f(0.5, -0.5, -0.5);      // P1 is red
	glColor3f(0.0, 1.0, 0.0);     glVertex3f(0.5, 0.5, -0.5);      // P2 is green
	glColor3f(0.0, 0.0, 1.0);     glVertex3f(-0.5, 0.5, -0.5);      // P3 is blue
	glColor3f(1.0, 0.0, 1.0);     glVertex3f(-0.5, -0.5, -0.5);      // P4 is purple

	glEnd();
	
	// Multi-colored plane
	glBegin(GL_POLYGON);
		glColor3f(1.0, 1.0, 1.0);
		glVertex3f(0.5, -0.5, 0.5);
		glVertex3f(-0.5, -0.5, 0.5);
		glVertex3f(-0.5, -0.5, -0.5);
		glVertex3f(0.5, -0.5, -0.5);
	glEnd();

	/*float dt = 0.01f;  //simulation time
	float tini = 0.0f; 
	float tfinal = 6.0f; //final time of simulation 
	// particle inicialization
	Particle p(0.0f, 5.0f, 0.0f); //initial position of the particle

	p.setLifetime(20.0f);
	std::cout << "Lifetime =" << p.getLifetime() << std::endl;
	p.setBouncing(0.8f);
	p.addForce(0, -9.8f, 0);
//	p.setFixed(true);

	// define a floor plane for collision
	glm::vec3 punt(0.0f);
	glm::vec3 normal(0, 1, 0);
	Plane plane(punt, normal);
	// simulation loop
	int count = 0;
	float disact, disant;
	disact = plane.distPoint2Plane(p.getCurrentPosition());
	float time = tini;

	while (time <= tfinal)
	{

		if (p.getLifetime() > 0) {
			p.updateParticle(dt, Particle::UpdateMethod::EulerOrig);
			std::cout << "position = " << p.getCurrentPosition().x << "  " << p.getCurrentPosition().y
				<< "  " << p.getCurrentPosition().z << "  time = " << time << std::endl;
			//Check for floor collisions
			disant = disact;
			p.setLifetime(p.getLifetime() - dt);
			disact = plane.distPoint2Plane(p.getCurrentPosition());
			if (disant*disact < 0.0f) { 
				//VERY IMPORTANT: only valid for the plane y=0 (floor plane)
				//Must be addapted to a general plane,
				p.setPosition(p.getCurrentPosition().x, -p.getCurrentPosition().y, p.getCurrentPosition().z);
				p.setVelocity(p.getVelocity().x, -p.getVelocity().y, p.getVelocity().z);
				std::cout << "Bounce = " << count++ << std::endl;
				disact = -disact; //
				system("PAUSE");
			}
		}
		time = time + dt; //increase time counter
	}
	system("PAUSE");
	return;*/

	glFlush();
	glutSwapBuffers();
}
