#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <glm\glm.hpp>
#include <vector>
#include "cal3d\cal3d.h"

class Particle
{
public:
	enum class UpdateMethod : std::int8_t { EulerOrig, EulerSemi, Verlet };

	Particle();
	Particle(const float& x, const float& y, const float& z);
//	Particle(glm::vec3 pos, glm::vec3 vel, float bouncing = 1.0f, bool fixed = false, int lifetime = -1, glm::vec3 force = glm::vec3(0, 0, 0));
	~Particle();
	//setters
	void setPosition(const float& x, const float& y, const float& z);
	void setPosition(glm::vec3 pos);
	void setPreviousPosition(const float& x, const float& y, const float& z);
	void setPreviousPosition(glm::vec3 pos);
	void setVelocity(const float& x, const float& y, const float& z);
	void setVelocity(glm::vec3 vel);
	void setForce(const float& x, const float& y, const float& z);
	void setForce(glm::vec3 force);
	void setBouncing(float bouncing);
	void setLifetime(float lifetime);
	void setFixed(bool fixed);
	void setWaypoint(glm::vec2 wp);
	void setNextWaypoint();

	//getters
	glm::vec3 getCurrentPosition();
	glm::vec3 getPreviousPosition();
	glm::vec3 getForce();
	glm::vec3 getVelocity();
	float getBouncing();
	float getLifetime();
	bool isFixed();
	float getDistanceToWaypoint();
	glm::vec3 getWaypoint();
	

	//other
	void addForce(glm::vec3 force);
	void addForce(const float& x, const float& y, const float& z);
	void updateParticle(const float& dt, UpdateMethod method = UpdateMethod::EulerOrig);
	void addStringfForce(glm::vec3 force);
	void clearStringForces();

	bool aStar(glm::vec2 start, glm::vec2 goal);

	int inOpen(glm::vec2 id);
	int inClose(glm::vec2 id);

	std::vector<glm::vec3> road;
	float max_speed;

	glm::vec3 avoidance;

private:
	glm::vec3 m_currentPosition;
	glm::vec3 m_previousPosition;
	glm::vec3 m_force;
	glm::vec3 m_velocity;
	glm::vec3 desired_velocity;
	std::vector<glm::vec3> stringForces;

	CalCoreModel* m_calCoreModel;
	CalModel* m_calModel;

	std::vector<std::vector<float>> gScores;
	std::vector<std::vector<float>> hScores;
	std::vector<std::vector<glm::vec2>> path;
	std::vector<glm::vec2> open;
	std::vector<glm::vec2> close;
	std::vector<std::vector<bool>> map;

	glm::vec3 waypoint;

	float m_bouncing;
	float m_lifetime;
	bool  m_fixed;

	float moveCost(glm::vec2 current, glm::vec2 n);
	std::vector<glm::vec2> getNeighboursOf(glm::vec2 idx);
	glm::vec3 getLowest();

};

#endif