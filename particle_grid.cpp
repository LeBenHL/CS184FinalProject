#include "particle_grid.h"
#include <iostream>
#include <algorithm>

ParticleGrid::ParticleGrid(ThreeDVector* min_bounds, ThreeDVector* max_bounds) {
    extern float H;
    this->min_bounds = min_bounds;

    int x = ceil((max_bounds->x - min_bounds->x)/H);
    int y = ceil((max_bounds->y - min_bounds->y)/H);
    int z = ceil((max_bounds->z - min_bounds->z)/H);

    this->grid_size = new ThreeDVector(x, y, z);
    this->neighbors_map = new map<ThreeDVector*, vector<Particle*>*, comparator>;
    this->water_particles = new vector<Particle*>;
    this->fog_particles = new vector<Particle*>;
    this->boundary_particles = new vector<Particle*>;

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
    delete this->water_particles;
    delete this->fog_particles;
    delete this->boundary_particles;
}

bool ParticleGrid::withInGrid(int x, int y, int z) {
    return x >= 0 && x < this->grid_size->x && y >= 0 && y < this->grid_size->y && z >= 0 && z < this->grid_size->z;
}

void ParticleGrid::addToGrid(Particle* particle) {
    switch (particle->type) {
        case Particle_Water: {
            this->water_particles->push_back(particle);
            break;
        }
        case Particle_Fog: {
            this->fog_particles->push_back(particle);
            break;
        }
        case Particle_Boundary: {
            this->boundary_particles->push_back(particle);
            break;
        }
    }
    this->registerGridPos(particle);
}

bool ParticleGrid::registerGridPos(Particle* particle) {
    extern float H;
    int x = (particle->position->x - this->min_bounds->x) / H;  
    int y = (particle->position->y - this->min_bounds->y) / H; 
    int z = (particle->position->z - this->min_bounds->z) / H; 

    if (withInGrid(x, y, z)) {
        #pragma omp critical(addToGrid) 
        (this->grid)[x][y][z]->push_back(particle);
        return true;
    } else {
        //Remove in iterator?
        this->removeFromGrid(particle, false);
        return false;
    }
}

bool ParticleGrid::unregisterGridPos(Particle* particle) {
    extern float H;
    int x = (particle->position->x - this->min_bounds->x) / H;  
    int y = (particle->position->y - this->min_bounds->y) / H; 
    int z = (particle->position->z - this->min_bounds->z) / H;

    if (withInGrid(x, y, z)) {
        vector<Particle*>* particles_cell =  (this->grid)[x][y][z];
        #pragma omp critical(removeFromGrid) 
        particles_cell->erase(std::remove(particles_cell->begin(), particles_cell->end(), particle), particles_cell->end());
        return true;
    } else {
        return false;
    }
}

void ParticleGrid::removeFromGrid(Particle* particle, bool unregister) {
    if (unregister) {
        this->unregisterGridPos(particle);
    }
    switch (particle->type) {
        case Particle_Water: {
            #pragma omp critical(eraseWater) 
            this->water_particles->erase(std::remove(this->water_particles->begin(), this->water_particles->end(), particle), this->water_particles->end());
            break;
        }
        case Particle_Fog: {
            #pragma omp critical(eraseFog) 
            this->fog_particles->erase(std::remove(this->fog_particles->begin(), this->fog_particles->end(), particle), this->fog_particles->end());
            break;
        }
        case Particle_Boundary: {
            #pragma omp critical(eraseBoundary) 
            this->boundary_particles->erase(std::remove(this->boundary_particles->begin(), this->boundary_particles->end(), particle), this->boundary_particles->end());
            break;
        }
    }
}

vector<Particle*>* ParticleGrid::getNeighbors(Particle* particle) {
    return getNeighbors(particle->position->x, particle->position->y, particle->position->z);
}

vector<Particle*>* ParticleGrid::getNeighbors(float pos_x, float pos_y, float pos_z) {
    extern float H;
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
        #pragma omp critical(neighborsMapInsert) 
        this->neighbors_map->insert(pair<ThreeDVector*, vector<Particle*>*>(position, neighbors));
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
    extern float H;
    return new ThreeDVector(x * H + this->min_bounds->x, y * H + this->min_bounds->y, z * H + this->min_bounds->z);
}