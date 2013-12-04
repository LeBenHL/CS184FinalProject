
#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <vector>
#include "three_d_vector.h"

using namespace std;

class Particle {
	public:
		ThreeDVector* position;
		long double mass;
		ThreeDVector* velocity;
		ThreeDVector* velocity_half;
		long double viscosity_coefficient;
		long double gas_constant;
		long double buoyancy_strength;
		long double rest_density;
		long double temperature;

		//Interpolated Fields
		long double density;
		ThreeDVector* acceleration;

		Particle(long double x, long double y, long double z, long double mass, 
			ThreeDVector* velocity, long double viscosity_coefficient, long double buoyancy_strength, 
			long double gas_constant, long double rest_density, long double temperature);
		~Particle();
		void set_density(vector<Particle*> particles);
		void set_acceleration(vector<Particle*> particles);

		//Leapfrog steps from http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf
		void leapfrog_start(long double dt);
		void leapfrog_step(long double dt);
		
		
		ThreeDVector* viscosityForce(vector<Particle*> particles);
		ThreeDVector* pressureForce(vector<Particle*> particles);
		ThreeDVector* pressure();
		ThreeDVector* externalForce();
		ThreeDVector* gravity();
		ThreeDVector* wind();
		ThreeDVector* buoyancy();

		long double color(vector<Particle*> particles);

		static Particle* createWaterParticle(long double x, long double y, long double z, ThreeDVector* velocity=new ThreeDVector());
		static Particle* createFogParticle(long double x, long double y, long double z, ThreeDVector* velocity=new ThreeDVector());


		//Poly 6 Kernel from http://www.matthiasmueller.info/publications/sca03.pdf
		static long double poly6Kernel(long double r, long double h);

		//Vicosity Kernel from http://www.matthiasmueller.info/publications/sca03.pdf
		static long double viscosityGradientSquaredKernel(long double r, long double h);

		//Spiky Kernel from http://www.matthiasmueller.info/publications/sca03.pdf
		static long double spikyGradientKernel(long double r, long double h);

};

#endif