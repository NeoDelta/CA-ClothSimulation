#include <iostream>
#include "Particle.h"
#include <glm\glm.hpp>
#include <vector>

class Constraint
{
	private:
		float rest_distance; // the length between particle p1 and p2 in rest configuration
		float elasticity, damping;

	public:
		Particle *p1, *p2; // the two particles that are connected through this constraint

		Constraint(Particle& par1, Particle& par2, float elast, float damp);

		/* This is one of the important methods, where a single constraint between two particles p1 and p2 is solved
		the method is called by Cloth.time_step() many times per frame*/
		void satisfyConstraint();
		void calculateStringForces();

};