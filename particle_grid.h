#ifndef PARTICLEGRID_H
#define PARTICLEGRID_H

#include "particle.h"
#include "three_d_vector.h"

using namespace std;

class ParticleGrid {
	//Code From by http://paulbourke.net/geometry/polygonise/
	public:
		ThreeDVector* min_bounds;
		ThreeDVector* grid_size;
		vector<vector<vector<vector<Particle*>* > > > grid;

		ParticleGrid(ThreeDVector* min_bounds, ThreeDVector* max_bounds);
		~ParticleGrid();

		void addToGrid(Particle* particle);
		bool withInGrid(int x, int y, int z);
		vector<Particle*>* getNeighbors(Particle* particle);
		vector<Particle*>* getNeighbors(long double pos_x, long double pos_y, long double pos_z);
		void removeFromGrid(Particle* particle);

};

#endif