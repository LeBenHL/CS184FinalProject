#include <cmath>

#include "particle.h"

using namespace std;

Particle::Particle(long double x, long double y, long double z, long double mass, 
		long double velocity, long double viscosity_coefficient, long double buoyancy_strength, long double gas_constant, long double rest_density) {
	this->position = new ThreeDVector(x, y, z);
	this->mass = mass;
	this->velocity = velocity;
	this->viscosity_coefficient = viscosity_coefficient;
	this->buoyancy_strength = buoyancy_strength;
	this->gas_constant = gas_constant;
	this->rest_density = rest_density;
}

Particle::~Particle() {
	delete this->position;
}

long double Particle::density(vector<Particle*> particles) {
	long double running_sum = 0;
	for (vector<Particle*>::iterator it = particles.begin(); it != particles.end(); it++) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		running_sum += particle->mass * Particle::poly6Kernel(this->position->distance(particle->position), 1);
	}
	//TODO Save density?
	return running_sum;
}

long double Particle::viscosityForce(vector<Particle*> particles) {
	long double running_sum = 0;
	for (vector<Particle*>::iterator it = particles.begin(); it != particles.end(); it++) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		running_sum += (particle->velocity - this->velocity) * (particle->mass / particle->density(particles)) * Particle::viscosityGradientSquaredKernel(this->position->distance(particle->position), 1);
	}
	return running_sum * this->viscosity_coefficient;
}

long double Particle::pressureForce(vector<Particle*> particles) {
	long double running_sum = 0;
	for (vector<Particle*>::iterator it = particles.begin(); it != particles.end(); it++) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		running_sum += ((particle->pressure(particles) + this->pressure(particles)) / 2) * (particle->mass / particle->density(particles)) * Particle::spikyGradientKernel(this->position->distance(particle->position), 1);
	}
	return -running_sum;
}

long double Particle::pressure(vector<Particle*> particles) {
	return this->gas_constant * (pow(this->density(particles) / this->rest_density, 7) - 1);
}

Particle* Particle::createWaterParticle(long double x, long double y, long double z, long double velocity) {
	extern long double WATER_MASS;
	extern long double WATER_VICOSITY_COEFFICIENT;
	extern long double WATER_BUOYANCY_STRENGTH;
	extern long double WATER_GAS_CONSTANT;
	extern long double WATER_REST_DENSITY;
	return new Particle(x, y, z, WATER_MASS, velocity, WATER_VICOSITY_COEFFICIENT, WATER_BUOYANCY_STRENGTH, WATER_GAS_CONSTANT, WATER_REST_DENSITY);
}

Particle* Particle::createFogParticle(long double x, long double y, long double z, long double velocity) {
	extern long double FOG_MASS;
	extern long double FOG_VICOSITY_COEFFICIENT;
	extern long double FOG_BUOYANCY_STRENGTH;
	extern long double FOG_GAS_CONSTANT;
	extern long double FOG_REST_DENSITY;
	return new Particle(x, y, z, FOG_MASS, velocity, FOG_VICOSITY_COEFFICIENT, FOG_BUOYANCY_STRENGTH, FOG_GAS_CONSTANT, FOG_REST_DENSITY);
}

//Using a Gaussian Kernal Function
long double Particle::poly6Kernel(long double r, long double h) {
	if (r >= 0 && r <= h) {
		extern long double PI;
		return (315 * pow(pow(h, 2) - pow(r, 2), 3)) / (64 * PI * pow(h, 9));
	} else {
		return 0;
	}
}

long double Particle::viscosityGradientSquaredKernel(long double r, long double h) {
	if (r >= 0 && r <= h) {
		extern long double PI;
		return (45 * (h - r)) / (PI * pow(h, 6));
	} else {
		return 0;
	}
}

long double Particle::spikyGradientKernel(long double r, long double h) {
	if (r >= 0 && r <= h) {
		extern long double PI;
		return -((45 * pow(h - r, 2)) / (PI * pow(h, 6)));
	} else {
		return 0;
	}
}