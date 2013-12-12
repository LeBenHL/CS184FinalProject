
#ifndef MARCHINGCUBE_H
#define MARCHINGCUBE_H

#include <vector>
#include "three_d_vector.h"
#include "particle.h"
#include "particle_grid.h"

using namespace std;

class MarchingCube {
	//Code From by http://paulbourke.net/geometry/polygonise/
	public:
		ThreeDVector* min_corner;
		float size;

		MarchingCube(ThreeDVector* min_corner, float size);
		~MarchingCube();
		vector<vector<pair<ThreeDVector*, ThreeDVector*> > >* triangulate(ParticleGrid* particle_grid, float isovalue_threshold);
		pair<ThreeDVector*, ThreeDVector*> interpolatePoint(int p1_index, int p2_index, float color_p1, float color_p2, float isovalue_threshold, ParticleGrid* particle_grid);
		ThreeDVector* pointAt(int corner_index);

		static vector<MarchingCube*>* generateGrid(vector<Particle*>* particles, float step_size);
		//static vector<MarchingCube*>* generateGridFast(ParticleGrid* particle_grid, float step_size);

		static const int edgeTable[256];
		static const int triTable[256][16];

};

#endif