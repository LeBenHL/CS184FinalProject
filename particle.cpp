#include <cmath>
#include "particle.h"
#include <iostream>
//#include <unordered_map>

using namespace std;

//CONSTANTS
float p_1 = 16769023.0;
float p_2 = 83492791.0;
float p_3 = 73856093.0;

//hash function

//prime multiplication method
/*unsigned long int hash_simple(ThreeDVector point, ) {
  float x = point->x * p_1;
  float y = point->y * p_2;
  float z = point->z * p_3;

  return (unsigned long int) (x ^ y ^ z);
}
*/

struct Key {
  ThreeDVector* point;

  bool operator==(const Key &other) const {
    return (point->x == other.point->x && point->y == other.point->y && point->z == other.point->z);
  }

};

/*
namespace std {
  template <>
  struct hash<Key> {
  public:
    std::size_t operator()(Key const* k) const {
      using std::size_t;
      using std::hash;

      float x = k->point->x * p_1;
      float y = k->point->y * p_2;
      float z = k->point->z * p_3;

      std::size_t h1 = std::hash<float>()(x);
      std::size_t h2 = std::hash<float>()(y);
      std::size_t h3 = std::hash<float>()(z);

      return h1 ^ h2 ^ h3; 
    }
  };
}*/

Particle::Particle(float x, float y, float z, float mass, 
		ThreeDVector* velocity, float viscosity_coefficient, float buoyancy_strength, 
		float gas_constant, float rest_density, float temperature, Particle_Type t, ThreeDVector* normal) {
	this->position = new ThreeDVector(x, y, z);
	this->mass = mass;
	this->velocity = velocity;
	this->viscosity_coefficient = viscosity_coefficient;
	this->buoyancy_strength = buoyancy_strength;
	this->gas_constant = gas_constant;
	this->rest_density = rest_density;
	this->temperature = temperature;
	this->acceleration = new ThreeDVector(0, 0, 0);
	this->type = t;
	this->normal = normal;
}

Particle::~Particle() {
	delete this->position;
	delete this->velocity;
	delete this->velocity_half;
	delete this->acceleration;
	delete this->normal;
}

void Particle::set_density(vector<Particle*>* particles) {
	extern float H;
	float running_sum = 0;
	if(this->type != Particle_Boundary){
		for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
			Particle* particle = *it;
			if(this->type == particle->type){
				running_sum += particle->mass * Particle::poly6Kernel(this->position->distance(particle->position), H);
			}
		}
	}
	//TODO Save density?

	//cout << running_sum << endl;
	this->density = running_sum;
}

void Particle::set_acceleration(vector<Particle*>* particles) {
	//Assumes we have previously called set_density for all particles to have density values
	//to use for this leapfrog step

	if(this->type != Particle_Boundary){
		ThreeDVector* net_force = new ThreeDVector();

		this->addExternalForce(net_force);
		this->addPressureForce(particles, net_force);
		this->addViscosityForce(particles, net_force);
		this->applyBoundaryConditions(particles, net_force);

		// Acceleration = Force / density
		net_force->scalar_multiply_bang(1 / this->density);

		delete this->acceleration;
		this->acceleration = net_force;
	}
}

void Particle::leapfrog_start(float dt) {
	if(this->type != Particle_Boundary){
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
}

void Particle::leapfrog_step(float dt) {
	if(this->type != Particle_Boundary){
		//Assumes that we have already calculated acceleration for this timestep
		//This is used to set half step velocities after the initial leapfrog

		//Set Half Step Velocity
		ThreeDVector a_dt = ThreeDVector(this->acceleration->x * dt, this->acceleration->y * dt, this->acceleration->z * dt);
		this->velocity_half->vector_add_bang(&a_dt);

		//Set Full Step Velocity For acceleration calculations
		ThreeDVector a_dt_over_2 = ThreeDVector(this->acceleration->x * dt/2, this->acceleration->y * dt/2, this->acceleration->z * dt/2);
		delete this->velocity;
		this->velocity = this->velocity_half->vector_add(&a_dt_over_2);

		//Set Full Step Position
		ThreeDVector v_dt = ThreeDVector(this->velocity_half->x * dt, this->velocity_half->y * dt, this->velocity_half->z * dt);
		this->position->vector_add_bang(&v_dt);
	}
}

void Particle::addViscosityForce(vector<Particle*>* particles, ThreeDVector* net_force) {
	extern float H;
	ThreeDVector running_sum = ThreeDVector();
	if (this->type != Particle_Boundary){
		for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
			Particle* particle = *it;
			if(this->type == particle->type){
				//TODO what is correct value of H?
				float multiplier = (particle->mass / particle->density) * Particle::viscosityGradientSquaredKernel(this->position->distance(particle->position), H);
				if (multiplier != 0) {
					ThreeDVector velocity_difference = ThreeDVector(particle->velocity->x - this->velocity->x, particle->velocity->y - this->velocity->y, particle->velocity->z - this->velocity->z);
					velocity_difference.scalar_multiply_bang(multiplier);
					running_sum.vector_add_bang(&velocity_difference);
				}
			}
		}
		running_sum.scalar_multiply_bang(this->viscosity_coefficient);
	}
	net_force->vector_add_bang(&running_sum);
}

void Particle::addPressureForce(vector<Particle*>* particles, ThreeDVector* net_force) {
	extern float H;
	ThreeDVector running_sum = ThreeDVector();
	if (this->type != Particle_Boundary) {
		float my_pressure = this->pressure();
		for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
			Particle* particle = *it;
			if(this->type == particle->type){
				//TODO what is correct value of H?
				ThreeDVector gradient;
				if (Particle::spikyGradientKernel(this->position, particle->position, H, &gradient)) {
					float particle_pressure = particle->pressure();
					float average_pressure = (particle_pressure + my_pressure) * 0.5 * (particle->mass / particle->density);
					gradient.scalar_multiply_bang(average_pressure);
					running_sum.vector_add_bang(&gradient);
				}
			}
		}
		running_sum.scalar_multiply_bang(-1);
	}
	net_force->vector_add_bang(&running_sum);
}

float parametric_calculation(float q){
	if (q > 0.0 && q < 2.0/3.0){
		return 2.0/3.0;
	}else if (q > 2.0/3.0 && q < 1.0){
		return 2.0*q - (3.0/2.0)*q*q;
	}else if(q > 1.0 && q < 2.0){
		return (1.0/2.0)*(2.0-q)*(2.0-q);
	}else{
		return 0.0;
	}
}

void Particle::addBoundaryForce(vector<Particle*>* particles, ThreeDVector* net_force) {
	extern float H;
	ThreeDVector running_sum = ThreeDVector();
	if(this->type != Particle_Boundary){
		for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
			Particle* particle = *it;
			if(particle->type == Particle_Boundary){
				float mass_a = this->mass;
				float mass_k = particle->mass;

				ThreeDVector xa_minus_xk = ThreeDVector(this->position->x - particle->position->x, this->position->y - particle->position->y, this->position->z - particle->position->z);
				float mag_xa_minus_xk = xa_minus_xk.magnitude();

				float q = mag_xa_minus_xk / H;
				
				float c2 = 100;

				float gamma = (0.02*c2)/mag_xa_minus_xk;
				gamma *= parametric_calculation(q);

				float constant_factor = ((mass_k/(mass_a + mass_k)) * gamma * this->density) / mag_xa_minus_xk;

				xa_minus_xk.scalar_multiply_bang(constant_factor);
				running_sum.vector_add_bang(&xa_minus_xk);
				//cout << c2 << endl;
				//cout << (mass_k/(mass_a + mass_k)) << endl;
				//cout << gamma << endl;
				//cout << mag_xa_minus_xk << endl;
				//cout << xa_minus_xk->repr() << endl;
				//cout << f_force->repr() << endl; 

			}
		}
	}
	net_force->vector_add_bang(&running_sum);
}

void Particle::applyBoundaryConditions(vector<Particle*>* particles, ThreeDVector* net_force) {
	float threshold = 0.1;
	if(this->type != Particle_Boundary){
		for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
			Particle* particle = *it;
			if(particle->type == Particle_Boundary && this->position->distance(particle->position) < threshold) {
				net_force->subtract_normal_component_bang(particle->normal);
				float multiplier = 100 * this->density;
				ThreeDVector boundary_force = ThreeDVector(particle->normal->x * multiplier, particle->normal->y * multiplier, particle->normal->z * multiplier);
				net_force->vector_add_bang(&boundary_force);
				this->velocity->subtract_normal_component_bang(particle->normal);
				this->velocity_half->subtract_normal_component_bang(particle->normal);
			}
		}
	}
}

void Particle::addExternalForce(ThreeDVector* net_force) {
	ThreeDVector external_force = ThreeDVector();
	this->addGravity(&external_force);
	this->addWind(&external_force);
	this->addBuoyancy(&external_force);

	net_force->vector_add_bang(&external_force);
}

void Particle::addGravity(ThreeDVector* vector) {
	extern ThreeDVector* CONSTANT_OF_GRAVITY;
	ThreeDVector gravity = ThreeDVector(CONSTANT_OF_GRAVITY->x * this->density, CONSTANT_OF_GRAVITY->y * this->density, CONSTANT_OF_GRAVITY->z * this->density);
	vector->vector_add_bang(&gravity);
}

void Particle::addWind(ThreeDVector* vector) {
	ThreeDVector wind = ThreeDVector();
	vector->vector_add_bang(&wind);
}

float Particle::pressure() {
	//Assume Pressure is same in all directions?
	float pressure = this->gas_constant * (pow(this->density / this->rest_density, 7) - 1);
	//cout << this->density / this->rest_density << endl;
	return pressure;
}

float ambientTemp(float height) {
	extern float AMBIENT_TEMP_AT_GROUND_LEVEL;
	return AMBIENT_TEMP_AT_GROUND_LEVEL + height * .3;//0.6
}

void Particle::addBuoyancy(ThreeDVector* vector) {
	extern ThreeDVector* CONSTANT_OF_GRAVITY;
	float multiplier = -this->buoyancy_strength * (this->temperature - ambientTemp(this->position->y)) * this->density;
	ThreeDVector buoyancy = ThreeDVector(CONSTANT_OF_GRAVITY->x * multiplier, CONSTANT_OF_GRAVITY->y * multiplier, CONSTANT_OF_GRAVITY->z * multiplier);
	vector->vector_add_bang(&buoyancy);
}

/*
bool Particle::isSurfaceParticle(vector<Particle*>* particles) {
	float delta = 0.05;
	float color = this->colorGradient(particles);
	cout << color << endl;
	return (color > (0.5 - delta)) && (color < (0.5 + delta));
}*/

float Particle::color(vector<Particle*>* particles) {
	return Particle::colorAt(this->position, particles, this->type);
}

/*
float Particle::colorGradient(vector<Particle*>* particles) {
	float step_size = 0.1;
	float x = this->position->x;
	float y = this->position->y;
	float z = this->position->z;

	float graident_x = (this->colorAt(x + step_size, y, z, particles) - this->colorAt(x - step_size, y, z, particles)) / step_size;
	float graident_y = (this->colorAt(x, y + step_size, z, particles) - this->colorAt(x - step_size, y - step_size, z, particles)) / step_size;
	float graident_z = (this->colorAt(x, y, z + step_size, particles) - this->colorAt(x - step_size, y, z - step_size, particles)) / step_size;
	ThreeDVector* gradient_vec = new ThreeDVector(graident_x, graident_y, graident_z);
	float magnitude = gradient_vec->magnitude();
	delete gradient_vec;
	return magnitude;
}*/

float Particle::colorAt(float x, float y, float z, vector<Particle*>* particles, Particle_Type t) {
	ThreeDVector position = ThreeDVector(x, y, z);
	float color = Particle::colorAt(&position, particles, t);
	return color;
}

float Particle::colorAt(ThreeDVector* position, vector<Particle*>* particles, Particle_Type t) {
	extern float H;
	map<ThreeDVector*, float, comparator>* color_map;
	if (t == Particle_Water) {
		color_map = Particle::water_color_map;
	} else {
		color_map = Particle::fog_color_map;
	}
	map<ThreeDVector*, float>::const_iterator got = color_map->find(position);
	if ( got == color_map->end() ) {
		float running_sum = 0;
		for (vector<Particle*>::iterator it = particles->begin(); it != particles->end(); ++it) {
			Particle* particle = *it;
			if (particle->type == t) {
				running_sum += particle->mass / particle->density * Particle::poly6Kernel(position->distance(particle->position), H);
			}
		}
		#pragma omp critical(colorMapInsert) 
		color_map->insert(pair<ThreeDVector*, float>(position->clone(), running_sum));
		return running_sum;
	} else {
		return got->second;
	}
}

Particle* Particle::createWaterParticle(float x, float y, float z, ThreeDVector* velocity) {
	extern float WATER_MASS;
	extern float WATER_VICOSITY_COEFFICIENT;
	extern float WATER_BUOYANCY_STRENGTH;
	extern float WATER_GAS_CONSTANT;
	extern float WATER_REST_DENSITY;
	extern float WATER_TEMP;
	return new Particle(x, y, z, WATER_MASS, velocity, WATER_VICOSITY_COEFFICIENT, WATER_BUOYANCY_STRENGTH, WATER_GAS_CONSTANT, WATER_REST_DENSITY, WATER_TEMP, Particle_Water);
}

Particle* Particle::createFogParticle(float x, float y, float z, ThreeDVector* velocity) {
	extern float FOG_MASS;
	extern float FOG_VICOSITY_COEFFICIENT;
	extern float FOG_BUOYANCY_STRENGTH;
	extern float FOG_GAS_CONSTANT;
	extern float FOG_REST_DENSITY;
	extern float FOG_TEMP;
	return new Particle(x, y, z, FOG_MASS, velocity, FOG_VICOSITY_COEFFICIENT, FOG_BUOYANCY_STRENGTH, FOG_GAS_CONSTANT, FOG_REST_DENSITY, FOG_TEMP, Particle_Fog);
}

Particle* Particle::createBoundaryParticle(float x, float y, float z, ThreeDVector* normal, ThreeDVector* velocity) {
	extern float BOUNDARY_MASS;
	return new Particle(x, y, z, BOUNDARY_MASS, velocity, 0, 0, 0, 0, 0, Particle_Boundary, normal);
}

float Particle::poly6Kernel(float r, float h) {
	if (r >= 0 && r <= h) {
		extern float PI;
		float hsquare_minus_rsquare = h * h - r * r;
		return (315 * hsquare_minus_rsquare * hsquare_minus_rsquare * hsquare_minus_rsquare) / (64 * PI * pow(h, 9));
	} else {
		return 0;
	}
}

float Particle::viscosityGradientSquaredKernel(float r, float h) {
	if (r >= 0 && r <= h) {
		extern float PI;
		return (45 * (h - r)) / (PI * pow(h, 6));
	} else {
		return 0;
	}
}

bool Particle::spikyGradientKernel(ThreeDVector* r, ThreeDVector* r_particle, float h, ThreeDVector* gradient) {
	float mag = r->distance(r_particle);
	if (mag >= 0 && mag <= h) {
		ThreeDVector delta = ThreeDVector(r->x - r_particle->x, r->y - r_particle->y, r->z - r_particle->z);
		extern float PI;
		//Normalize Vector
		if (mag != 0) {
			delta.scalar_multiply_bang(1/mag);
		}
		float h_minus_mag = h - mag;
		delta.scalar_multiply_bang(-((45 * h_minus_mag * h_minus_mag) / (PI * pow(h, 6))));
		gradient->x = delta.x;
		gradient->y = delta.y;
		gradient->z = delta.z;
		return true;
	} else {
		return false;
	}
}

void Particle::clearColorMap() {
	for(map<ThreeDVector*, float, comparator>::iterator it = Particle::water_color_map->begin(); it != Particle::water_color_map->end(); it++) {
	    delete it->first;
	}
	for(map<ThreeDVector*, float, comparator>::iterator it = Particle::fog_color_map->begin(); it != Particle::fog_color_map->end(); it++) {
	    delete it->first;
	}
	Particle::water_color_map->clear();
	Particle::fog_color_map->clear();
}

map<ThreeDVector*, float, comparator>* Particle::water_color_map = new map<ThreeDVector*, float, comparator>;
map<ThreeDVector*, float, comparator>* Particle::fog_color_map = new map<ThreeDVector*, float, comparator>;







