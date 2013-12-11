#include <cmath>

#include "particle.h"
#include <iostream>

using namespace std;

Particle::Particle(long double x, long double y, long double z, long double mass, 
		ThreeDVector* velocity, long double viscosity_coefficient, long double buoyancy_strength, 
		long double gas_constant, long double rest_density, long double temperature) {
	this->position = new ThreeDVector(x, y, z);
	this->mass = mass;
	this->velocity = velocity;
	this->viscosity_coefficient = viscosity_coefficient;
	this->buoyancy_strength = buoyancy_strength;
	this->gas_constant = gas_constant;
	this->rest_density = rest_density;
	this->temperature = temperature;
	this->acceleration = new ThreeDVector(0, 0, 0);
}

Particle::~Particle() {
	delete this->position;
	delete this->velocity;
	delete this->velocity_half;
	delete this->acceleration;
}

void Particle::set_density(vector<Particle*>* particles) {
	extern long double H;
	long double running_sum = 0;
	for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
		Particle* particle = *it;
		running_sum += particle->mass * Particle::poly6Kernel(this->position->distance(particle->position), H);
	}
	//TODO Save density?

	//cout << running_sum << endl;
	this->density = running_sum;
}

void Particle::set_acceleration(vector<Particle*>* particles) {
	//Assumes we have previously called set_density for all particles to have density values
	//to use for this leapfrog step

	ThreeDVector* net_force = new ThreeDVector();

	ThreeDVector* external_force = this->externalForce();
	ThreeDVector* pressure_force = this->pressureForce(particles);
	ThreeDVector* viscosity_force = this->viscosityForce(particles);

	net_force->vector_add_bang(external_force);
	net_force->vector_add_bang(pressure_force);
	net_force->vector_add_bang(viscosity_force);

	delete external_force;
	delete pressure_force;
	delete viscosity_force;

	// Acceleration = Force / density
	net_force->scalar_multiply_bang(1 / this->density);

	delete this->acceleration;
	this->acceleration = net_force;
}

void Particle::leapfrog_start(long double dt) {
	//Assumes that we have already calculated acceleration for this timestep
	//This is used to initialize the half time step velocities

	//Set Half Step Velocity
	ThreeDVector* a_dt_over_2 = this->acceleration->scalar_multiply(dt/2);
	this->velocity_half = this->velocity->vector_add(a_dt_over_2);
	delete a_dt_over_2;

	//Set Full Step Velocity
	ThreeDVector* a_dt = this->acceleration->scalar_multiply(dt);
	this->velocity->vector_add_bang(a_dt);
	delete a_dt;

	//Set Full Step Position
	ThreeDVector* v_dt = this->velocity_half->scalar_multiply(dt);
	this->position->vector_add_bang(v_dt);
	delete v_dt;
}

void Particle::leapfrog_step(long double dt) {
	//Assumes that we have already calculated acceleration for this timestep
	//This is used to set half step velocities after the initial leapfrog

	//Set Half Step Velocity
	ThreeDVector* a_dt = this->acceleration->scalar_multiply(dt);
	this->velocity_half->vector_add_bang(a_dt);
	delete a_dt;

	//Set Full Step Velocity For acceleration calculations
	ThreeDVector* a_dt_over_2 = this->acceleration->scalar_multiply(dt/2);
	delete this->velocity;
	this->velocity = this->velocity_half->vector_add(a_dt_over_2);
	delete a_dt_over_2;

	//Set Full Step Position
	ThreeDVector* v_dt = this->velocity_half->scalar_multiply(dt);
	this->position->vector_add_bang(v_dt);
	delete v_dt;
}

ThreeDVector* Particle::viscosityForce(vector<Particle*>* particles) {
	extern long double H;
	ThreeDVector* running_sum = new ThreeDVector();
	for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		ThreeDVector* velocity_difference = particle->velocity->vector_subtract(this->velocity);
		velocity_difference->scalar_multiply_bang(particle->mass / particle->density);
		velocity_difference->scalar_multiply_bang(Particle::viscosityGradientSquaredKernel(this->position->distance(particle->position), H));
		running_sum->vector_add_bang(velocity_difference);
		delete velocity_difference;
	}
	running_sum->scalar_multiply_bang(this->viscosity_coefficient);
	return running_sum;
}

ThreeDVector* Particle::pressureForce(vector<Particle*>* particles) {
	extern long double H;
	ThreeDVector* running_sum = new ThreeDVector();
	long double my_pressure = this->pressure();
	for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
		Particle* particle = *it;
		//TODO what is correct value of H?
		long double particle_pressure = particle->pressure();
		long double average_pressure = (particle_pressure + average_pressure) * 0.5 * (particle->mass / particle->density);
		ThreeDVector* r_hat = this->position->vector_subtract(particle->position);
		ThreeDVector* gradient = Particle::spikyGradientKernel(r_hat, H);
		gradient->scalar_multiply_bang(average_pressure);
		running_sum->vector_add_bang(gradient);
		delete r_hat;
		delete gradient;
	}

	running_sum->scalar_multiply_bang(-1);
	return running_sum;
}

ThreeDVector* Particle::externalForce() {
	ThreeDVector* external_force = new ThreeDVector();
	ThreeDVector* gravity = this->gravity();
	ThreeDVector* wind = this->wind();
	ThreeDVector* buoyancy = this->buoyancy();

	external_force->vector_add_bang(gravity);
	external_force->vector_add_bang(wind);
	external_force->vector_add_bang(buoyancy);

	delete gravity;
	delete wind;
	delete buoyancy;

	return external_force;
}

ThreeDVector* Particle::gravity() {
	extern ThreeDVector* CONSTANT_OF_GRAVITY;
	return CONSTANT_OF_GRAVITY->scalar_multiply(this->density);
}

ThreeDVector* Particle::wind() {
	return new ThreeDVector(5, 0, 5);
}

long double Particle::pressure() {
	//Assume Pressure is same in all directions?
	long double pressure = this->gas_constant * (pow(this->density / this->rest_density, 7) - 1);
	return pressure;
}

ThreeDVector* Particle::buoyancy() {
	extern long double AMBIENT_TEMP;
	extern ThreeDVector* CONSTANT_OF_GRAVITY;
	return CONSTANT_OF_GRAVITY->scalar_multiply(-this->buoyancy_strength * (this->temperature - AMBIENT_TEMP));
}

/*
bool Particle::isSurfaceParticle(vector<Particle*>* particles) {
	long double delta = 0.05;
	long double color = this->colorGradient(particles);
	cout << color << endl;
	return (color > (0.5 - delta)) && (color < (0.5 + delta));
}*/

long double Particle::color(vector<Particle*>* particles) {
	return Particle::colorAt(this->position, particles);
}

/*
long double Particle::colorGradient(vector<Particle*>* particles) {
	long double step_size = 0.1;
	long double x = this->position->x;
	long double y = this->position->y;
	long double z = this->position->z;

	long double graident_x = (this->colorAt(x + step_size, y, z, particles) - this->colorAt(x - step_size, y, z, particles)) / step_size;
	long double graident_y = (this->colorAt(x, y + step_size, z, particles) - this->colorAt(x - step_size, y - step_size, z, particles)) / step_size;
	long double graident_z = (this->colorAt(x, y, z + step_size, particles) - this->colorAt(x - step_size, y, z - step_size, particles)) / step_size;
	ThreeDVector* gradient_vec = new ThreeDVector(graident_x, graident_y, graident_z);
	long double magnitude = gradient_vec->magnitude();
	delete gradient_vec;
	return magnitude;
}*/

long double Particle::colorAt(long double x, long double y, long double z, vector<Particle*>* particles) {
	ThreeDVector* position = new ThreeDVector(x, y, z);
	long double color = Particle::colorAt(position, particles);
	delete position;
	return color;
}

long double Particle::colorAt(ThreeDVector* position, vector<Particle*>* particles) {
	extern long double H;
	map<ThreeDVector*, long double>::const_iterator got = Particle::color_map->find(position);
	if ( got == Particle::color_map->end() ) {
		long double running_sum = 0;
		for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
			Particle* particle = *it;
			running_sum += particle->mass / particle->density * Particle::poly6Kernel(position->distance(particle->position), H);
		}
		Particle::color_map->insert(pair<ThreeDVector*, long double>(position->clone(), running_sum));
		return running_sum;
	} else {
		return got->second;
	}
}

Particle* Particle::createWaterParticle(long double x, long double y, long double z, ThreeDVector* velocity) {
	extern long double WATER_MASS;
	extern long double WATER_VICOSITY_COEFFICIENT;
	extern long double WATER_BUOYANCY_STRENGTH;
	extern long double WATER_GAS_CONSTANT;
	extern long double WATER_REST_DENSITY;
	extern long double WATER_TEMP;
	return new Particle(x, y, z, WATER_MASS, velocity, WATER_VICOSITY_COEFFICIENT, WATER_BUOYANCY_STRENGTH, WATER_GAS_CONSTANT, WATER_REST_DENSITY, WATER_TEMP);
}

Particle* Particle::createFogParticle(long double x, long double y, long double z, ThreeDVector* velocity) {
	extern long double FOG_MASS;
	extern long double FOG_VICOSITY_COEFFICIENT;
	extern long double FOG_BUOYANCY_STRENGTH;
	extern long double FOG_GAS_CONSTANT;
	extern long double FOG_REST_DENSITY;
	extern long double FOG_TEMP;
	return new Particle(x, y, z, FOG_MASS, velocity, FOG_VICOSITY_COEFFICIENT, FOG_BUOYANCY_STRENGTH, FOG_GAS_CONSTANT, FOG_REST_DENSITY, FOG_TEMP);
}

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

ThreeDVector* Particle::spikyGradientKernel(ThreeDVector* r_hat, long double h) {
	long double r = r_hat->magnitude();
	if (r >= 0 && r <= h) {
		extern long double PI;
		return r_hat->scalar_multiply(-((45 * pow(h - r, 2)) / (PI * pow(h, 6))));
	} else {
		return new ThreeDVector();
	}
}

void Particle::clearColorMap() {
	for(map<ThreeDVector*, long double, comparator>::iterator it = Particle::color_map->begin(); it != Particle::color_map->end(); it++) {
	    delete it->first;
	}
	Particle::color_map->clear();
}

map<ThreeDVector*, long double, comparator>* Particle::color_map = new map<ThreeDVector*, long double, comparator>;






