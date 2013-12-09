
#ifndef MARCHINGCUBE_H
#define MARCHINGCUBE_H

#include <vector>
#include "three_d_vector.h"
#include "particle.h"

using namespace std;

class MarchingCube {
	//Code From by http://paulbourke.net/geometry/polygonise/
	public:
		ThreeDVector* min_corner;
		long double size;

		MarchingCube(ThreeDVector* min_corner, long double size);
		~MarchingCube();
		vector<vector<pair<ThreeDVector*, ThreeDVector*> > >* triangulate(vector<Particle*>* particles, long double isovalue_threshold);
		pair<ThreeDVector*, ThreeDVector*> interpolatePoint(int p1_index, int p2_index, long double color_p1, long double color_p2, long double isovalue_threshold, vector<Particle*>* particles);
		ThreeDVector* pointAt(int corner_index);

		static vector<MarchingCube*>* generateGrid(vector<Particle*>* particles, long double step_size);

		static const int edgeTable[256];
		static const int triTable[256][16];

};

#endif