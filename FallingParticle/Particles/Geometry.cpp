#pragma once
#include "Geometry.h"

//****************************************************
// Plane
//****************************************************

Plane::Plane() {
	normal = glm::normalize(glm::vec3(0.0, 1.0, 0.0));
	dconst = -glm::dot(glm::vec3(0.0, 0.0, 0.0), normal);
};

Plane::~Plane() {
};

Plane::Plane(const glm::vec3& point, const glm::vec3& normalVect){
	normal = glm::normalize(normalVect);
	dconst = -glm::dot(point, normal);
};

Plane::Plane(const glm::vec3& point0, const glm::vec3& point1, const glm::vec3& point2){
	glm::vec3 v1 = point1 - point0;
	glm::vec3 v2 = point2 - point0;
	normal = glm::normalize(glm::cross(v1, v2));
	dconst = -glm::dot(point0, normal);
};

void Plane::setPosition(const glm::vec3& newPos){
	dconst = -glm::dot(newPos, normal);
};

bool Plane::isInside(const glm::vec3& point){
	float dist;
	dist = glm::dot(point, normal) + dconst;
	if (dist > 1.e-7)
		return false;
	else
		return true;
};

float Plane::distPoint2Plane(const glm::vec3& point){
	float dist;
	return dist = glm::dot(point, normal) + dconst;
};

glm::vec3 Plane::closestPointInPlane(const glm::vec3& point){
	glm::vec3 closestP;
	float r = (-dconst - glm::dot(point, normal));
	return closestP = point + r*normal;
};

bool Plane::intersecSegment(const glm::vec3& point1, const glm::vec3& point2, glm::vec3& pTall){
	if (distPoint2Plane(point1)*distPoint2Plane(point2) > 0)	return false;
	float r = (-dconst - glm::dot(point1, normal)) / glm::dot((point2 - point1), normal);
	pTall = (1 - r)*point1 + r*point2;
	return true;
};


//****************************************************
// Triangle
//****************************************************

Triangle::Triangle(const glm::vec3& point0, const glm::vec3& point1, const glm::vec3& point2) {
	vertex1 = point0;
	vertex2 = point1;
	vertex3 = point2;

	normal = glm::normalize(glm::cross(vertex2 - vertex1, vertex3 - vertex1));
}


void Triangle::setPosition(const glm::vec3& newPos) {
	vertex1 = newPos;
	vertex2 = newPos + (vertex2 - vertex1);
	vertex3 = newPos + (vertex3 - vertex1);
}

bool Triangle::isInside(const glm::vec3& point) {
	return false;
}

bool Triangle::intersecSegment(const glm::vec3& point1, const glm::vec3& point2, glm::vec3& pTall) {

	glm::vec3 Normal, IntersectPos;


	// Find distance from LP1 and LP2 to the plane defined by the triangle
	float Dist1 = dot(point1 - vertex1, normal);
	float Dist2 = dot(point2 - vertex1, normal);

	if ((Dist1 * Dist2) >= 0.0f) return false; // line doesn't cross the triangle.

	if (Dist1 == Dist2)  return false; // line and plane are parallel

	  // Find point on the line that intersects with the plane
	IntersectPos = point1 + (vertex2 - point1) * (-Dist1 / (Dist2 - Dist1));

	// Find if the interesection point lies inside the triangle by testing it against all edges
	glm::vec3 vTest;

	vTest = cross(normal, vertex2 - vertex1);
	if (dot(vTest, IntersectPos - vertex1) < 0.0f) return false;

	vTest = cross(normal, vertex3 - vertex2);
	if (dot(vTest, IntersectPos - vertex2) < 0.0f) return false;

	vTest = cross(normal, vertex1 - vertex3);
	if (dot(vTest, IntersectPos - vertex1) < 0.0f) return false;

	pTall = IntersectPos;

	return true;
}

//****************************************************
// Sphere
//****************************************************

Sphere::Sphere(const glm::vec3& point, const float& radious) {
	center = point;
	radi = radious;
}

Sphere::~Sphere() {

}

void Sphere::setPosition(const glm::vec3& newPos) {
	center = newPos;
}

bool Sphere::isInside(const glm::vec3& point) {
	if (distPointCenter(point) <= radi) return true;
	return false;
}

float Sphere::distPointCenter(const glm::vec3& point) {
	glm::vec3 distVector = point - center;
	float dist = glm::length(distVector); //distVector.length();

	return dist;
}
