#include "particle_grid.h"
#include <iostream>

ParticleGrid::ParticleGrid(ThreeDVector* min_bounds, ThreeDVector* max_bounds) {
	extern long double H;
	this->min_bounds = min_bounds;

	int x = ceil((max_bounds->x - min_bounds->x)/H);
	int y = ceil((max_bounds->y - min_bounds->y)/H);
	int z = ceil((max_bounds->z - min_bounds->z)/H);

	this->grid_size = new ThreeDVector(x, y, z);
    this->neighbors_map = new map<ThreeDVector*, vector<Particle*>*, comparator>;

	this->grid.resize(x);

	for (int i = 0; i < x; i++) {
		this->grid[i].resize(y);
		for (int j = 0; j < y; j++) {
		  	this->grid[i][j].resize(z);
		  	for (int k = 0; k < z; k++) {
		    	this->grid[i][j][k] = new vector<Particle*>;
		  	}
		}
	}
}

ParticleGrid::~ParticleGrid() {
	delete this->min_bounds;
	delete this->grid_size;
    delete this->neighbors_map;
}

bool ParticleGrid::withInGrid(int x, int y, int z) {
    return x >= 0 && x < this->grid_size->x && y >= 0 && y < this->grid_size->y && z >= 0 && z < this->grid_size->z;
}

void ParticleGrid::addToGrid(Particle* particle) {
	extern long double H;
	int x = (particle->position->x - this->min_bounds->x) / H;  
	int y = (particle->position->y - this->min_bounds->y) / H; 
	int z = (particle->position->z - this->min_bounds->z) / H; 

	if (withInGrid(x, y, z)) {
	   (this->grid)[x][y][z]->push_back(particle);
	} else {
	}
}

void ParticleGrid::removeFromGrid(Particle* particle) {
	extern long double H;
	int x = (particle->position->x - this->min_bounds->x) / H;  
	int y = (particle->position->y - this->min_bounds->y) / H; 
	int z = (particle->position->z - this->min_bounds->z) / H;

	if (withInGrid(x, y, z)) {
    	vector<Particle*>* particles_cell =  (this->grid)[x][y][z];
    	particles_cell->erase(remove(particles_cell->begin(), particles_cell->end(), particle), particles_cell->end());
	}
}

vector<Particle*>* ParticleGrid::getNeighbors(Particle* particle) {
    return getNeighbors(particle->position->x, particle->position->y, particle->position->z);
}

vector<Particle*>* ParticleGrid::getNeighbors(long double pos_x, long double pos_y, long double pos_z) {
    extern long double H;
    int x = (pos_x - this->min_bounds->x) / H;  
    int y = (pos_y - this->min_bounds->y) / H; 
    int z = (pos_z - this->min_bounds->z) / H;

    ThreeDVector* position = new ThreeDVector(x, y, z);
    map<ThreeDVector*, vector<Particle*>*>::const_iterator got = this->neighbors_map->find(position);
    if ( got == this->neighbors_map->end() ) { 
        vector<Particle*>* neighbors = new vector<Particle*>;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                for (int k = -1; k <= 1; k++) {
                    if (withInGrid(x+i, y+j, z+k)) {
                      (*neighbors).insert((*neighbors).end(), (this->grid)[x+i][y+j][z+k]->begin(), (this->grid)[x+i][y+j][z+k]->end());
                    }
                }
            }
        }
        this->neighbors_map->insert(make_pair<ThreeDVector*, vector<Particle*>*>(position, neighbors));
        return neighbors;
    } else {
        delete position;
        return got->second;
    }
}

void ParticleGrid::clearNeighborsMap() {
    for(map<ThreeDVector*, vector<Particle*>*>::iterator it = this->neighbors_map->begin(); it != this->neighbors_map->end(); it++) {
        delete it->first;
        delete it->second;
    }
    this->neighbors_map->clear();
}

ThreeDVector* ParticleGrid::minCornerOfCell(int x, int y, int z) {
	extern long double H;
	return new ThreeDVector(x * H + this->min_bounds->x, y * H + this->min_bounds->y, z * H + this->min_bounds->z);
}