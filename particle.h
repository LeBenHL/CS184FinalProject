
#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <vector>
#include <map>
#include <iostream>
#include "three_d_vector.h"

using namespace std;

typedef enum {
  Particle_Water,
  Particle_Fog,
  Particle_Boundary
}Particle_Type;

class Particle {
	public:
		ThreeDVector* position;
		float mass;
		ThreeDVector* velocity;
		ThreeDVector* velocity_half;
		float viscosity_coefficient;
		float gas_constant;
		float buoyancy_strength;
		float rest_density;
		float temperature;

		 Particle_Type type;

		static map<ThreeDVector*, float, comparator>* color_map;

		//Interpolated Fields
		float density;
		ThreeDVector* acceleration;

		Particle(float x, float y, float z, float mass, 
			ThreeDVector* velocity, float viscosity_coefficient, float buoyancy_strength, 
			float gas_constant, float rest_density, float temperature, Particle_Type t);
		~Particle();
		void set_density(vector<Particle*>* particles);
		void set_acceleration(vector<Particle*>* particles);

		//Leapfrog steps from http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf
		void leapfrog_start(float dt);
		void leapfrog_step(float dt);
		
		
		ThreeDVector* viscosityForce(vector<Particle*>* particles);
		ThreeDVector* pressureForce(vector<Particle*>* particles);
		ThreeDVector* boundaryForce(vector<Particle*>* particles);
		float pressure();
		ThreeDVector* externalForce();
		void addGravity(ThreeDVector* vector);
		void addWind(ThreeDVector* vector);
		void addBuoyancy(ThreeDVector* vector);

		//bool isSurfaceParticle(vector<Particle*>* particles);

		float color(vector<Particle*>* particles);
		//float colorGradient(vector<Particle*>* particles);

		static float colorAt(ThreeDVector* position, vector<Particle*>* particles);
		static float colorAt(float x, float y, float z, vector<Particle*>* particles);

		static Particle* createWaterParticle(float x, float y, float z, ThreeDVector* velocity=new ThreeDVector());
		static Particle* createFogParticle(float x, float y, float z, ThreeDVector* velocity=new ThreeDVector());
		static Particle* createBoundaryParticle(float x, float y, float z, ThreeDVector* velocity=new ThreeDVector());


		//Poly 6 Kernel from http://www.matthiasmueller.info/publications/sca03.pdf
		static float poly6Kernel(float r, float h);

		//Vicosity Kernel from http://www.matthiasmueller.info/publications/sca03.pdf
		static float viscosityGradientSquaredKernel(float r, float h);

		//Spiky Kernel from http://www.matthiasmueller.info/publications/sca03.pdf
		static bool spikyGradientKernel(ThreeDVector* r, ThreeDVector* r_particle, float h, ThreeDVector* gradient);

		static void clearColorMap();

};

#endif