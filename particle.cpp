#include "particle.h"

using namespace std;

Particle::Particle(long double x, long double y, long double z, long double mass, 
		long double velocity, long double viscosity_coefficient, long double buoyancy_strength) {
	this->position = new ThreeDVector(x, y, z);
	this->mass = mass;
	this->velocity = velocity;
	this->viscosity_coefficient = viscosity_coefficient;
	this->buoyancy_strength = buoyancy_strength;
}

Particle::~Particle() {
	delete this->position;
}

Particle* Particle::createWaterParticle(long double x, long double y, long double z, long double velocity) {
	long double WATER_MASS = 1.0;
	long double WATER_VICOSITY_COEFFICIENT = 1.0;
	long double WATER_BUOYANCY_STRENGTH = 1.0;
	return new Particle(x, y, z, WATER_MASS, velocity, WATER_VICOSITY_COEFFICIENT, WATER_BUOYANCY_STRENGTH);
}

Particle* Particle::createFogParticle(long double x, long double y, long double z, long double velocity) {
	long double FOG_MASS = 1.0;
	long double FOG_VICOSITY_COEFFICIENT = 1.0;
	long double FOG_BUOYANCY_STRENGTH = 1.0;
	return new Particle(x, y, z, FOG_MASS, velocity, FOG_VICOSITY_COEFFICIENT, FOG_BUOYANCY_STRENGTH);
}