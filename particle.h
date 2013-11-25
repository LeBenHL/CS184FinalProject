#ifndef PARTICLE_H
#define PARTICLE_H

#include "three_d_vector.h"

class Particle {
	public:
		ThreeDVector* position;
		long double mass;
		long double velocity;
		long double viscosity_coefficient;
		long double buoyancy_strength;

		Particle(long double x, long double y, long double z, long double mass, 
			long double velocity, long double viscosity_coefficient, long double buoyancy_strength);
		~Particle();
		static Particle* createWaterParticle(long double x, long double y, long double z, long double velocity=0);
		static Particle* createFogParticle(long double x, long double y, long double z, long double velocity=0);

};

#endif