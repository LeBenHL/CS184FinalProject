#include <cmath>

#include "particle.h"

using namespace std;

Particle::Particle(long double x, long double y, long double z, long double mass, 
		ThreeDVector* velocity, long double viscosity_coefficient, long double buoyancy_strength, long double gas_constant, long double rest_density) {
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
	delete this->velocity;
}

void Particle::set_density(vector<Particle*> particles) {
	long double running_sum = 0;
	for (vector<Particle*>::iterator it = particles.begin(); it != particles.end(); it++) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		running_sum += particle->mass * Particle::poly6Kernel(this->position->distance(particle->position), 1);
	}
	//TODO Save density?
	this->density = running_sum;
}

ThreeDVector* Particle::viscosityForce(vector<Particle*> particles) {
	ThreeDVector* running_sum = new ThreeDVector();
	for (vector<Particle*>::iterator it = particles.begin(); it != particles.end(); it++) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		ThreeDVector* velocity_difference = particle->velocity->vector_subtract(this->velocity);
		velocity_difference->scalar_multiply_bang(particle->mass / particle->density);
		velocity_difference->scalar_multiply_bang(Particle::viscosityGradientSquaredKernel(this->position->distance(particle->position), 1));
		running_sum->vector_add_bang(velocity_difference);
		delete velocity_difference;
	}
	running_sum->scalar_multiply_bang(this->viscosity_coefficient);
	return running_sum;
}

ThreeDVector* Particle::pressureForce(vector<Particle*> particles) {
	ThreeDVector* running_sum = new ThreeDVector();
	ThreeDVector* my_pressure = this->pressure();
	for (vector<Particle*>::iterator it = particles.begin(); it != particles.end(); it++) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		ThreeDVector* particle_pressure = particle->pressure();
		ThreeDVector* average_pressure = particle_pressure->vector_add(my_pressure);
		average_pressure->scalar_multiply_bang(0.5);
		average_pressure->scalar_multiply_bang(particle->mass / particle->density);
		average_pressure->scalar_multiply_bang(Particle::spikyGradientKernel(this->position->distance(particle->position), 1));
		running_sum->vector_add_bang(average_pressure);
		delete particle_pressure;
		delete average_pressure;
	}
	delete my_pressure;

	running_sum->scalar_multiply_bang(-1);
	return running_sum;
}

ThreeDVector* Particle::externalForce() {
	ThreeDVector* externalForce = new ThreeDVector();
	ThreeDVector* gravity = this->gravity();
	ThreeDVector* wind = this->wind();

	externalForce->vector_add_bang(gravity);
	externalForce->vector_add_bang(wind);

	delete gravity;
	delete wind;

	return externalForce;
}

ThreeDVector* Particle::gravity() {
	extern ThreeDVector* CONSTANT_OF_GRAVITY;
	return CONSTANT_OF_GRAVITY->scalar_multiply(this->density);
}

ThreeDVector* Particle::wind() {
	return new ThreeDVector(5, 0, 0);
}

ThreeDVector* Particle::pressure() {
	//Assume Pressure is same in all directions?
	long double pressure = this->gas_constant * (pow(this->density / this->rest_density, 7) - 1);
	return new ThreeDVector(pressure, pressure, pressure);
}

Particle* Particle::createWaterParticle(long double x, long double y, long double z, ThreeDVector* velocity) {
	extern long double WATER_MASS;
	extern long double WATER_VICOSITY_COEFFICIENT;
	extern long double WATER_BUOYANCY_STRENGTH;
	extern long double WATER_GAS_CONSTANT;
	extern long double WATER_REST_DENSITY;
	return new Particle(x, y, z, WATER_MASS, velocity, WATER_VICOSITY_COEFFICIENT, WATER_BUOYANCY_STRENGTH, WATER_GAS_CONSTANT, WATER_REST_DENSITY);
}

Particle* Particle::createFogParticle(long double x, long double y, long double z, ThreeDVector* velocity) {
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