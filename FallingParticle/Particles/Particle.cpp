#if defined(_MSC_VER) && _MSC_VER <= 0x0600
#pragma warning(disable : 4786)
#endif

#include "Particle.h"


Particle::Particle()
{
	glm::vec3 pos(0, 0, 0);
	m_currentPosition = pos;
	m_fixed = false;

	//m_calCoreModel = new CalCoreModel("dummy");
}

Particle::Particle(const float& x, const float& y, const float& z) :
m_previousPosition(0, 0, 0), m_velocity(0, 0, 0), m_force(0, 0, 0), m_bouncing(1), m_lifetime(50), m_fixed(false), desired_velocity(0, 0, 0), max_speed(0), avoidance(0, 0, 0)
{
	m_currentPosition.x = x;
	m_currentPosition.y = y;
	m_currentPosition.z = z;

	for (int i = 0; i < 20; i++) {
		std::vector<bool> aux;
		for (int j = 0; j < 20; j++) {
			if (i == 10 && j < 9) aux.push_back(false);
			else if (i == 10 && j > 11) aux.push_back(false);
			else aux.push_back(true);
		}
		map.push_back(aux);
	}

}

/*
Particle::Particle(glm::vec3 pos, glm::vec3 vel, float bouncing, bool fixed, int lifetime, glm::vec3 force) :
m_currentPosition(pos), m_previousPosition(pos), m_force(force), m_velocity(vel), m_bouncing(bouncing), m_lifetime(lifetime), m_fixed(fixed)
{
}
*/

Particle::~Particle()
{
}

//setters
void Particle::setPosition(const float& x, const float& y, const float& z)
{
	glm::vec3 pos(x,y,z);
	m_currentPosition =  pos;
}
void Particle::setPosition(glm::vec3 pos)
{
	m_currentPosition = pos;
}

void Particle::setPreviousPosition(const float& x, const float& y, const float& z)
{
	glm::vec3 pos(x, y, z);
	m_previousPosition = pos;
}

void Particle::setPreviousPosition(glm::vec3 pos)
{
	m_previousPosition = pos;
}

void Particle::setForce(const float& x, const float& y, const float& z)
{
	glm::vec3 force(x, y, z);
	m_force = force;
}

void Particle::setForce(glm::vec3 force)
{
	m_force = force;
}

void Particle::addForce(const float& x, const float& y, const float& z)
{
	glm::vec3 force(x,y,z);
	m_force += force;
}

void Particle::addForce(glm::vec3 force)
{
	m_force += force;
}

void Particle::setVelocity(const float& x, const float& y, const float& z)
{
	glm::vec3 vel(x,y,z);
	m_velocity = vel;
}

void Particle::setVelocity(glm::vec3 vel)
{
	m_velocity = vel;
}

void Particle::setBouncing(float bouncing)
{
	m_bouncing = bouncing;
}

void Particle::setLifetime(float lifetime)
{
	m_lifetime = lifetime;
}

void Particle::setFixed(bool fixed)
{
	m_fixed = fixed;
}

//getters
glm::vec3 Particle::getCurrentPosition()
{
	return m_currentPosition;
}

glm::vec3 Particle::getPreviousPosition()
{
	return m_previousPosition;
}

glm::vec3 Particle::getForce()
{
	return m_force;
}

glm::vec3 Particle::getVelocity()
{
	return m_velocity;
}

float Particle::getBouncing()
{
	return m_bouncing;
}

float Particle::getLifetime()
{
	return m_lifetime;
}

bool Particle::isFixed()
{
	return m_fixed;
}

void Particle::addStringfForce(glm::vec3 force)
{
	stringForces.push_back(force);
}

void Particle::clearStringForces()
{
	stringForces.clear();
}

void Particle::updateParticle(const float& dt, UpdateMethod method)
{
	glm::vec3 force = glm::vec3(0.0f);
	float damp = 0.01;
	if (!m_fixed & m_lifetime > 0)
	{
		for (int f = 0; f < stringForces.size(); ++f) {
			force += stringForces.at(f);
		}
		switch (method)
		{
		case UpdateMethod::EulerOrig:
		{
			m_previousPosition = m_currentPosition;
			m_currentPosition += m_velocity*dt;
			m_velocity += (m_force + force)*dt;
		}
			break;
		case UpdateMethod::EulerSemi:
		{
			glm::vec3 p1_to_p2 = waypoint - m_currentPosition;
			float norm = sqrt(p1_to_p2.x * p1_to_p2.x + p1_to_p2.y * p1_to_p2.y + p1_to_p2.z * p1_to_p2.z); // current distance between p1 and p2
			p1_to_p2 = p1_to_p2 / norm;

			desired_velocity = p1_to_p2 * max_speed;
			glm::vec3 steering = desired_velocity - m_velocity;
			
			if (avoidance != glm::vec3(0, 0, 0)) {
				m_velocity = avoidance * max_speed;
			}
			else {
				m_velocity += steering*15.0f*dt;
			}
			//m_velocity += steering;

			m_previousPosition = m_currentPosition;
			m_velocity += (m_force + force)*dt;
			m_currentPosition += m_velocity*dt;

			avoidance = glm::vec3(0, 0, 0);

		}
			break;
		case UpdateMethod::Verlet:
		{
			glm::vec3 aux = m_currentPosition;

			//m_currentPosition += (m_currentPosition - m_previousPosition)*(1 - damp) + (m_force + force)*dt*dt;
			m_currentPosition += (m_currentPosition - m_previousPosition) + (m_force + force)*dt*dt;
			m_velocity = (m_currentPosition - aux) / dt;

			m_previousPosition = aux;
		}
		break;
		}
	}

	if (!stringForces.empty()) stringForces.clear();
	return;
}

bool Particle::aStar(glm::vec2 start, glm::vec2 goal)
{
	open.clear();
	close.clear();
	gScores.clear();
	hScores.clear();
	path.clear();
	road.clear();

	for (int i = 0; i < 20; i++) {
		std::vector<float> aux;
		std::vector<glm::vec2> aux2;
		for (int j = 0; j < 20; j++) {
			 aux.push_back(10000000);
			 aux2.push_back(glm::vec2(-1, -1));
		}
		gScores.push_back(aux);
		hScores.push_back(aux);
		path.push_back(aux2);
	}
	gScores[start[0]][start[1]] = 0;
	hScores[start[0]][start[1]] = moveCost(start, goal);

	open.push_back(start);
	glm::vec3 res = getLowest();

	glm::vec2 lowest = glm::vec2(res[0], res[1]);
	int openIdx = res[2];

	while (lowest != goal) {
		res = getLowest();

		lowest = glm::vec2(res[0], res[1]);
		openIdx = res[2];

		glm::vec2 current = lowest;
		open.erase(open.begin() + openIdx);
		close.push_back(current);

		std::vector<glm::vec2> neigh = getNeighboursOf(current);

		for (unsigned int i = 0; i < neigh.size(); i++) {
			glm::vec2 n = neigh[i];

			float cost = gScores[current[0]][current[1]] + moveCost(current, n);

			int nOpen = inOpen(n);
			int nClosed = inClose(n);

			if (nOpen != -1 && cost < gScores[n[0]][n[1]]) open.erase(open.begin() + nOpen);
			if (nClosed != -1 && cost < gScores[n[0]][n[1]]) close.erase(close.begin() + nClosed);

			if (nOpen == -1 && nClosed == -1) {

				gScores[n[0]][n[1]] = cost;
				hScores[n[0]][n[1]] = cost + moveCost(n, goal);

				open.push_back(n);
				path[n[0]][n[1]] = current;
			}
		}

	}

	//Path generator
	glm::vec2 current;
	current = goal;
	road.push_back(glm::vec3(goal[0]-9.5f, 0.25f, goal[1]-9.5f));
	while (current != start) {
		current = path[current[0]][current[1]];
		road.push_back(glm::vec3(current[0]-9.5f, 0.25f, current[1]-9.5f));
	} 

	if (!road.empty()) 
		waypoint = road.back();
		road.pop_back();

	if (!road.empty()) {
		waypoint = road.back();

		/*glm::vec3 p1_to_p2 = waypoint - m_currentPosition;
		float norm = sqrt(p1_to_p2.x * p1_to_p2.x + p1_to_p2.y * p1_to_p2.y + p1_to_p2.z * p1_to_p2.z); // current distance between p1 and p2
		p1_to_p2 = p1_to_p2 / norm;

		m_velocity = p1_to_p2 * max_speed;*/
		road.pop_back();
	}
	

	return true;
}

glm::vec3 Particle::getLowest()
{
	float min = 100000.0f;
	glm::vec3 res;

	for (unsigned int i = 0; i < open.size(); i++) {
		glm::vec2 idx = open[i];
		if (hScores[idx[0]][idx[1]] < min) {
			min = hScores[idx[0]][idx[1]];
			res = glm::vec3(idx[0], idx[1], i);
		}
	}

	return res;
}

std::vector<glm::vec2> Particle::getNeighboursOf(glm::vec2 idx)
{
	std::vector<glm::vec2> neigh;
	int x = idx[0];
	int y = idx[1];

	if (x > 0) 
		if (map[x - 1][y])	neigh.push_back(glm::vec2(x - 1, y)); 
	if (x < 19) 
		if (map[x + 1][y])	neigh.push_back(glm::vec2(x + 1, y));
	if (y > 0) 
		if (map[x][y - 1])	neigh.push_back(glm::vec2(x, y - 1)); 
	if (y < 19) 
		if( map[x][y + 1])	neigh.push_back(glm::vec2(x, y + 1));
	if (x > 0 && y > 0 )	
		if( map[x - 1][y - 1] && map[x][y - 1] && map[x -1][y])	neigh.push_back(glm::vec2(x - 1, y - 1));
	if (x < 19 && y < 19)	
		if( map[x + 1][y + 1] && map[x][y + 1] && map[x + 1][y])	neigh.push_back(glm::vec2(x + 1, y + 1));
	if (x < 19 && y > 0)	
		if( map[x + 1][y - 1] && map[x][y - 1] && map[x + 1][y])	neigh.push_back(glm::vec2(x + 1, y - 1));
	if (x > 0 && y < 19)	
		if( map[x - 1][y + 1] && map[x][y + 1] && map[x - 1][y])	neigh.push_back(glm::vec2(x - 1, y + 1));

	return neigh;
}

float Particle::moveCost(glm::vec2 current, glm::vec2 n)
{
	float xc = current[0];
	float yc = current[1];

	float xn = n[0];
	float yn = n[1];

	return sqrt(pow(xn - xc, 2) + pow(yn - yc, 2));
}

int Particle::inOpen(glm::vec2 id)
{
	for (int i = 0; i < open.size(); i++) {
		if (open[i] == id) return i;
	}

	return -1;
}

int Particle::inClose(glm::vec2 id)
{
	for (int i = 0; i < close.size(); i++) {
		if (close[i] == id) return i;
	}

	return -1;
}

glm::vec3 Particle::getWaypoint() 
{
	return waypoint;
}

float Particle::getDistanceToWaypoint()
{
	float px = m_currentPosition[0];
	float py = m_currentPosition[2];

	float wx = waypoint[0];
	float wy = waypoint[2];

	return sqrt(pow(wx - px, 2) + pow(wy - py, 2));
}

void Particle::setNextWaypoint()
{
	if (!road.empty()) {
		waypoint = road.back();

		/*glm::vec3 p1_to_p2 = waypoint - m_currentPosition;
		float norm = sqrt(p1_to_p2.x * p1_to_p2.x + p1_to_p2.y * p1_to_p2.y + p1_to_p2.z * p1_to_p2.z); // current distance between p1 and p2
		p1_to_p2 = p1_to_p2 / norm;

		m_velocity = p1_to_p2 * max_speed;*/
		road.pop_back();
	}
}