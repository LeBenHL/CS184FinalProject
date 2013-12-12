
#include <vector>
#include <queue>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <limits>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>

#include "three_d_vector.h"
#include "lodepng.h"
#include "particle.h"
#include "marching_cube.h"
#include "particle_grid.h"

//CONSTANTS
float PI = atan(1)*4;
float E = 2.7182818284590452353;
ThreeDVector* CONSTANT_OF_GRAVITY = new ThreeDVector(0, -9.8, 0);
float AMBIENT_TEMP = 25;
float TIMESTEP_DURATION =  0.0125;
float PARTICLE_RADIUS = 0.04;
//float H = 0.01;
float H = .225;
float MARCHING_CUBE_STEP_SIZE = .1;
//float MARCHING_CUBE_STEP_SIZE = 2;
float ISOVALUE_THRESHOLD = 0.5;

float WATER_MASS = 1.0;
float WATER_VICOSITY_COEFFICIENT = 3.5;
float WATER_BUOYANCY_STRENGTH = 0.01;
float WATER_GAS_CONSTANT = 500;
float WATER_REST_DENSITY = 650;
float WATER_TEMP = 30.0;
float FOG_MASS = 1.0;

static int img_counter = 0;

//best values so far
/*
ThreeDVector* CONSTANT_OF_GRAVITY = new ThreeDVector(0, -9.8, 0);
float AMBIENT_TEMP = 25;
float TIMESTEP_DURATION =  0.01;
float PARTICLE_RADIUS = 0.1;
//float H = 0.01;
float H = .225;
float MARCHING_CUBE_STEP_SIZE = .1;
//float MARCHING_CUBE_STEP_SIZE = 2;
float ISOVALUE_THRESHOLD = 0.5;

float WATER_MASS = 1.0;
float WATER_VICOSITY_COEFFICIENT = 1;
float WATER_BUOYANCY_STRENGTH = 0.01;
float WATER_GAS_CONSTANT = 500;
float WATER_REST_DENSITY = 700;
float WATER_TEMP = 30.0;
*/

float FOG_VICOSITY_COEFFICIENT = 1.0;
float FOG_BUOYANCY_STRENGTH = 1.0;
float FOG_GAS_CONSTANT = 1.0;
float FOG_REST_DENSITY = 1.0;
float FOG_TEMP = 1.0;
float BOUNDARY_MASS = 20;

//COLORS
typedef enum {
    Color_Golden_Gate_Orange,
    Color_Ground_Brown,
    Water,
    Fog
}Color;

using namespace std;

//Our Polygons/Primatives
vector<vector<pair<ThreeDVector*, ThreeDVector*> > > polygons;

//Save Boolean
bool save = false;
//Filename
char* file_name;

//Types of Surface Reconstruction
bool spheres = true;
bool marching_cubes = true;

//How Many Timesteps we have advanced so far
float num_timesteps = 0;

//Max/Min x,y,z
float max_x = numeric_limits<float>::min();
float min_x = numeric_limits<float>::max();
float max_y = numeric_limits<float>::min();
float min_y = numeric_limits<float>::max();
float max_z = numeric_limits<float>::min();
float min_z = numeric_limits<float>::max();

//Bounding Box
ThreeDVector* min_bounds;
ThreeDVector* max_bounds;

//Grids containing Particles
ParticleGrid* particle_grid;

float scale_factor = 1.0/500000.0;

//Print Function for debugging
void print(string _string) {
  cout << _string << endl;
}

void print(float num) {
  cout << num << endl;
}

void print(int num) {
  cout << num << endl;
}

//Create PNG Function
void createPng(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
  //Encode the image
  unsigned error = lodepng::encode(filename, image, width, height);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

void parseObj(const char* filename) {
  cout << "Parsing Object File" << endl;
  
  vector<ThreeDVector*> vertices;
  vector<ThreeDVector*> vertices_normals;

  std::ifstream inpfile(filename);
  if(!inpfile.is_open()) {
    std::cout << "Unable to open file" << std::endl;
  } else {
    std::string line;
    //MatrixStack mst;
    
    while(inpfile.good()) {
      std::vector<std::string> splitline;
      std::string buf;

      std::getline(inpfile,line);
      std::stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }
      //Valid commands:
      //v x y z [w]
      else if(!splitline[0].compare("v")) {
        float x = atof(splitline[1].c_str());
        float y = atof(splitline[2].c_str());
        float z = atof(splitline[3].c_str());
        ThreeDVector* vertex = new ThreeDVector(x, y, z);
        vertices.push_back(vertex);
      }
      //vn x y z
      else if(!splitline[0].compare("vn")) {
        float x = atof(splitline[1].c_str());
        float y = atof(splitline[2].c_str());
        float z = atof(splitline[3].c_str());
        ThreeDVector* normal = new ThreeDVector(x, y, z);
        normal->normalize_bang();
        vertices_normals.push_back(normal);
      }
      //f v1 v2 v3 v4 ....
      //We assume they are all triangles defined with 3 vertices
      //Most files seem to be like this
      //1 indexed vertices so we need to add 1
      else if(!splitline[0].compare("f")) {
        vector<pair<ThreeDVector*, ThreeDVector*> > polygon;

        for(int i = 1; i < splitline.size(); i++){
          const char* v = splitline[i].c_str();
          char v_str[500];
          strncpy(v_str, v, sizeof(v_str));
          vector<int> vector_of_indices; 
          char *pch = strtok(v_str, "/");
          int last = 0;
          while (pch != NULL){
            last++;
            vector_of_indices.push_back(atoi(pch) - 1);
            pch = strtok(NULL,"/");
          }  
          ThreeDVector* vertex = vertices[vector_of_indices[0]];
          ThreeDVector* normal = vertices_normals[vector_of_indices[last-1]];

          pair<ThreeDVector*, ThreeDVector*> pair = std::make_pair(vertex, normal);
          polygon.push_back(pair);

          max_x = max(max_x, vertex->x);
          min_x = min(min_x, vertex->x);
          max_y = max(max_y, vertex->y);
          min_y = min(min_y, vertex->y);
          max_z = max(max_z, vertex->z);
          min_z = min(min_z, vertex->z);
        }

        polygons.push_back(polygon);
      }
    }
  }
  
  vertices.clear();
}


//old setBounds
/*
void setBounds() {
  float center_x = (max_x + min_x) / 2;
  float center_y = (max_y + min_y) / 2;
  float center_z = (max_z + min_z) / 2;

  float transformed_max_x = (max_x - center_x) * scale_factor;
  float transformed_min_x = (min_x - center_x) * scale_factor;
  float transformed_max_y = (max_y - center_y) * scale_factor;
  float transformed_min_y = (min_y - center_y) * scale_factor;
  float transformed_max_z = (max_z - center_z) * scale_factor;
  float transformed_min_z = (min_z - center_z) * scale_factor;

  min_bounds = new ThreeDVector(transformed_min_x * 2, transformed_min_y * 2, transformed_min_z * 10);
  max_bounds = new ThreeDVector(transformed_max_x * 2, transformed_max_y * 2, transformed_max_z * 10);

  particle_grid = new ParticleGrid(min_bounds, max_bounds);
}*/

void setBounds() {
  min_bounds = new ThreeDVector(-3.0, -3.0, -3.0);
  max_bounds = new ThreeDVector(5.0, 5.0, 5.0);
  
  particle_grid = new ParticleGrid(min_bounds, max_bounds);
}

void advanceOneTimestep() {
  //Calculate densities for this time step first
  #pragma omp parallel for
  for (int i = 0; i < particle_grid->water_particles->size(); ++i) {
    Particle* particle = particle_grid->water_particles->at(i);
    vector<Particle*>* neighbors;
    neighbors = particle_grid->getNeighbors(particle);
    particle->set_density(neighbors);
    //delete neighbors;
  }

  //Calculate accelerations next
  #pragma omp parallel for
  for (int i = 0; i < particle_grid->water_particles->size(); ++i) {
    Particle* particle = particle_grid->water_particles->at(i);
    vector<Particle*>* neighbors = particle_grid->getNeighbors(particle);
    particle->set_acceleration(neighbors);
    //delete neighbors;
  }

  //Then perform leapfrog!
  vector<Particle*> water_copy(*(particle_grid->water_particles));
  #pragma omp parallel for
  for (int i = 0; i < water_copy.size(); ++i) {
    Particle* particle = water_copy.at(i);

    particle_grid->unregisterGridPos(particle);
    //First Timestep
    if (num_timesteps == 0) {
      particle->leapfrog_start(TIMESTEP_DURATION);
    //Not First Timestep
    } else {
      particle->leapfrog_step(TIMESTEP_DURATION);
    }
    particle_grid->registerGridPos(particle);
  }

  //Calculate densities for this time step first
  #pragma omp parallel for
  for (int i = 0; i < particle_grid->fog_particles->size(); ++i) {
    Particle* particle = particle_grid->fog_particles->at(i);
    particle->set_density(particle_grid->fog_particles);
  }

  //Calculate accelerations next
  #pragma omp parallel for
  for (int i = 0; i < particle_grid->fog_particles->size(); ++i) {
    Particle* particle = particle_grid->fog_particles->at(i);
    particle->set_acceleration(particle_grid->fog_particles);
  }

  //Then perform leapfrog!
  for (vector<Particle*>::iterator it = particle_grid->fog_particles->begin(); it != particle_grid->fog_particles->end(); ++it) {
    Particle* particle = *it;
    
    //First Timestep
    if (num_timesteps == 0) {
      particle->leapfrog_start(TIMESTEP_DURATION);
    //Not First Timestep
    } else {
      particle->leapfrog_step(TIMESTEP_DURATION);
    }
  }

  //Also clear color map since that info is not valid anymore!
  Particle::clearColorMap();
  //Clear Neighbors Maps as well!
  particle_grid->clearNeighborsMap();

  ++num_timesteps;
}

//****************************************************
// Some Classes
//****************************************************

class Viewport;

class Viewport {
  public:
    int w, h; // width and height
};


//****************************************************
// Global Variables
//****************************************************
Viewport  viewport;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  //glOrtho((center_x - radius) * multiplier, (center_x + radius) * multiplier, (center_y - radius) * multiplier, (center_y + radius) * multiplier, 1.0, 1.0 + radius);
  gluPerspective(90, float(w)/float(h), 1.0, 11.0);
  glViewport (0,0,viewport.w,viewport.h);

}

//****************************************************
// Simple init function
//****************************************************
void initScene(){
  //glClearColor(0.52941f, 0.80784f, 0.98039f, 0.0f); // Clear to Sky Blue, fully transparent
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  // Enable lighting and the light we have set up
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glEnable(GL_LIGHT1);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_RESCALE_NORMAL);

  //Set lighting parameters
  GLfloat light_position0[] = {0.0, 1.0, 1.0, 0};
  GLfloat light_ambient0[] = {0, 0, 0, 1};
  GLfloat light_diffuse0[] = {2.0, 2.0, 2.0, 1};
  GLfloat light_specular0[] = {2.0, 2.0, 2.0, 1};

  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

  glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.0);
  glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0);
  glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.0);

  GLfloat light_position1[] = {0.0, 0.0, 1.0, 0};
  GLfloat light_ambient1[] = {0, 0, 0, 1};
  GLfloat light_diffuse1[] = {1.0, 1.0, 1.0, 1};
  GLfloat light_specular1[] = {1.0, 1.0, 1.0, 1};

  glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
  glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);

  glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 1.0);
  glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.0);
  glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.0);

  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  myReshape(viewport.w,viewport.h);
}

void setColor(int color) {
  switch (color) {
    case Color_Golden_Gate_Orange: {
      GLfloat ambient_color[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat diffuse_color[] = { 0.753, 0.211, 0.173, 1.0 };
      GLfloat specular_color[] = { 0.3, 0.3, 0.3, 1.0 };
      GLfloat shininess[] = { 50.0 };
      GLfloat emission[] = {0, 0, 0, 1};

      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
      break;
    }
    case Color_Ground_Brown: {
      GLfloat ambient_color[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat diffuse_color[] = { 0.41176, 0.27451, 0.16078, 1.0 };
      GLfloat specular_color[] = { 0, 0, 0, 1.0 };
      GLfloat shininess[] = { 50.0 };
      GLfloat emission[] = {0, 0, 0, 1};

      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
      break;
    }
    case Water: {
      GLfloat ambient_color[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat diffuse_color[] = { 0.14510, 0.42745, 0.48235, 1.0 };
      GLfloat specular_color[] = { 0.3, 0.3, 0.3, 1.0 };
      GLfloat shininess[] = { 50.0 };
      GLfloat emission[] = {0, 0, 0, 1};

      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
      break;
    }
    case Fog: {
      GLfloat ambient_color[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat diffuse_color[] = { 204.0/255.0, 207.0/255.0, 188.0/255.0, 1.0 };
      GLfloat specular_color[] = { 0.3, 0.3, 0.3, 1.0 };
      GLfloat shininess[] = { 50.0 };
      GLfloat emission[] = {0, 0, 0, 1};

      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
      break;
    }
  }
}

void marchingCubes(vector<Particle*>* particles) {
  int counter = 0;
  vector<MarchingCube*>* cubes = MarchingCube::generateGrid(particles, MARCHING_CUBE_STEP_SIZE);
  //vector<MarchingCube*>* cubes = MarchingCube::generateGridFast(particle_grid, MARCHING_CUBE_STEP_SIZE);
  vector<vector<pair<ThreeDVector*, ThreeDVector*> > >* triangles = new vector<vector<pair<ThreeDVector*, ThreeDVector*> > >;

  #pragma omp parallel for
  for(int i = 0; i < cubes->size(); ++i) {
    MarchingCube* cube = cubes->at(i);
    //char* buffer = new char[1000];
    //sprintf(buffer, "%d/%d", ++counter, cubes->size());
    //print(buffer);
    vector<vector<pair<ThreeDVector*, ThreeDVector*> > >* triangles_subset = cube->triangulate(particle_grid, ISOVALUE_THRESHOLD);
    #pragma omp critical
    triangles->insert(triangles->end(), triangles_subset->begin(), triangles_subset->end());
    delete triangles_subset;
    delete cube;
  }

  for(vector<vector<pair<ThreeDVector*, ThreeDVector*> > >::iterator it = triangles->begin(); it != triangles->end(); ++it) {
    vector<pair<ThreeDVector*, ThreeDVector*> > triangle = *it;
    glBegin(GL_POLYGON);                      // Draw A Polygon
    for(vector<pair<ThreeDVector*, ThreeDVector*> >::iterator i = triangle.begin(); i != triangle.end(); ++i) {
      pair<ThreeDVector*, ThreeDVector*> vertex_pair = *i;
      ThreeDVector* vertex = vertex_pair.first;
      ThreeDVector* normal = vertex_pair.second;
      glNormal3f(normal->x, normal->y, normal->z);
      glVertex3f(vertex->x, vertex->y, vertex->z);
      delete vertex;
      delete normal;
    }
  glEnd();
  }

  delete cubes;
}


//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);       // clear the color buffer

  glMatrixMode(GL_MODELVIEW);         // indicate we are specifying camera transformations
  glLoadIdentity();                   // make sure transformation is "zero'd"

  float center_x = (max_x + min_x) / 2;
  float center_y = (max_y + min_y) / 2;
  float center_z = (max_z + min_z) / 2;
  //print((min_y - center_y) * scale_factor);

  //gluLookAt(-5, 0, 2.5, -4, 0, 0, 0, 1, 0);
  gluLookAt(-4, 3, 4, 0, 0, -1, 0, 1, 0);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Enable shading
  glShadeModel(GL_SMOOTH);
  /*
  //BRIDGE
  setColor(Color_Golden_Gate_Orange);
  glPushMatrix();
  glScalef(scale_factor, scale_factor, scale_factor);
  glTranslatef(-center_x, -center_y, -center_z);
  // Start drawing
  for(vector<vector<pair<ThreeDVector*, ThreeDVector*> > >::iterator it = polygons.begin(); it != polygons.end(); ++it) {
    vector<pair<ThreeDVector*, ThreeDVector*> > polygon = *it;
    glBegin(GL_POLYGON);                      // Draw A Polygon
    for(vector<pair<ThreeDVector*, ThreeDVector*> >::iterator i = polygon.begin(); i != polygon.end(); ++i) {
      pair<ThreeDVector*, ThreeDVector*> vertex_pair = *i;
      ThreeDVector* vertex = vertex_pair.first;
      ThreeDVector* normal = vertex_pair.second;
      glNormal3f(normal->x, normal->y, normal->z);
      glVertex3f(vertex->x, vertex->y, vertex->z);
    }
    glEnd();
  }
  glPopMatrix();

  //BOUNDING SURFACES
  setColor(Color_Ground_Brown);
  glBegin(GL_POLYGON);  
  glNormal3f(0, 1, 0);
  glVertex3f(-10000, -35.0, -10000);
  glNormal3f(0, 1, 0);
  glVertex3f(-10000, -35.0, 10000);
  glNormal3f(0, 1, 0);
  glVertex3f(10000, -35.0, 10000);
  glNormal3f(0, 1, 0);
  glVertex3f(10000, -35.0, -10000);
  glEnd();
  //FLUIDS

  setColor(Color_Ground_Brown);
  glBegin(GL_POLYGON);  
  glNormal3f(0, 1, 0);
  glVertex3f(-10000, -35.0, -10000);
  glNormal3f(0, 1, 0);
  glVertex3f(-10000, -35.0, 10000);
  glNormal3f(0, 1, 0);
  glVertex3f(10000, -35.0, 10000);
  glNormal3f(0, 1, 0);
  glVertex3f(10000, -35.0, -10000);
  glEnd();  */

  setColor(Color_Ground_Brown); //bottom wall
  glBegin(GL_POLYGON);  
  glNormal3f(0, 1, 0);
  glVertex3f(-3, -1, -3);
  glNormal3f(0, 1, 0);
  glVertex3f(-3, -1, 3);
  glNormal3f(0, 1, 0);
  glVertex3f(3, -1, 3);
  glNormal3f(0, 1, 0);
  glVertex3f(3, -1, -3);
  glEnd();

  glBegin(GL_LINE_LOOP); //left wall 
  glNormal3f(1, 0, 0);
  glVertex3f(-3, -1, -3);
  glNormal3f(1, 0, 0);
  glVertex3f(-3, 1, -3);
  glNormal3f(1, 0, 0);
  glVertex3f(-3, 1, 3);
  glNormal3f(1, 0, 0);
  glVertex3f(-3, -1, 3);
  glEnd();

  glBegin(GL_LINE_LOOP); //right wall 
  glNormal3f(-1, 0, 0);
  glVertex3f(3, -1, -3);
  glNormal3f(-1, 0, 0);
  glVertex3f(3, 1, -3);
  glNormal3f(-1, 0, 0);
  glVertex3f(3, 1, 3);
  glNormal3f(-1, 0, 0);
  glVertex3f(3, -1, 3);
  glEnd();

  glBegin(GL_LINE_LOOP); //front wall 
  glNormal3f(0, 0, -1);
  glVertex3f(-3, -1, 3);
  glNormal3f(0, 0, -1);
  glVertex3f(-3, 1, 3);
  glNormal3f(0, 0, -1);
  glVertex3f(3, 1, 3);
  glNormal3f(0, 0, -1);
  glVertex3f(3, -1, 3);
  glEnd();

  glBegin(GL_LINE_LOOP); //back wall 
  glNormal3f(0, 0, 1);
  glVertex3f(-3, -1, -3);
  glNormal3f(0, 0, 1);
  glVertex3f(-3, 1, -3);
  glNormal3f(0, 0, 1);
  glVertex3f(3, 1, -3);
  glNormal3f(0, 0, 1);
  glVertex3f(3, -1, -3);
  glEnd();

  if (spheres) {
    setColor(Water);
    for (vector<Particle*>::iterator it = particle_grid->water_particles->begin(); it != particle_grid->water_particles->end(); it++) {
      Particle* particle = *it;
      glPushMatrix();
      glTranslatef(particle->position->x, particle->position->y, particle->position->z);
      glutSolidSphere(PARTICLE_RADIUS, 20, 20);
      glPopMatrix();
    }

    setColor(Fog);
    for (vector<Particle*>::iterator it = particle_grid->fog_particles->begin(); it != particle_grid->fog_particles->end(); it++) {
      Particle* particle = *it;
      glPushMatrix();
      glTranslatef(particle->position->x, particle->position->y, particle->position->z);
      glutSolidSphere(PARTICLE_RADIUS, 20, 20);
      glPopMatrix();
    }
    /*
    setColor(Color_Ground_Brown);
    for (vector<Particle*>::iterator it = particle_grid->boundary_particles->begin(); it != particle_grid->boundary_particles->end(); it++) {
      Particle* particle = *it;
      glPushMatrix();
      glTranslatef(particle->position->x, particle->position->y, particle->position->z);
      glutSolidSphere(PARTICLE_RADIUS, 100, 100);
      glPopMatrix();
    }*/
  } else if (marching_cubes) {
    setColor(Water);
    marchingCubes(particle_grid->water_particles);

    setColor(Fog);
    marchingCubes(particle_grid->fog_particles);
  }



  if (save) {
    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    vector<unsigned char> buf(w * h * 4);

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, &buf[0] );

    std::ostringstream os_str;
    os_str << file_name << img_counter << ".png";

    cout << os_str.str().c_str() << endl;

    const char* f_name = os_str.str().c_str();

    img_counter++;
    createPng(f_name, buf, w, h);
  }

  glFlush();
  glutSwapBuffers();          // swap buffers (we earlier set double buffer)
}



//****************************************************
// the usual stuff, nothing exciting here
//****************************************************

void myKeyboardFunc(unsigned char key, int x, int y){
}

void mySpecialKeyFunc(int key, int x, int y){

}

void myFrameMove() {
  //nothing here for now
  #ifdef _WIN32
    Sleep(10);                                   //give ~10ms back to OS (so as not to waste the CPU)
  #endif
    advanceOneTimestep();
    glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

int main(int argc, char *argv[]) {

  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-save") {
      if(i + 1 < argc) {
        save = true;
        file_name = argv[i + 1];
        i = i + 1;
      }
    } 
  }

  //Parse Polygons the Golden Gate
  //parseObj("Golden Gate Bridge.obj");
  setBounds();
  
  for (int x = -13; x <13; x++) {
    for (int y = -3; y <9; y++) {
      for (int z = -13; z <13; z++) {
        Particle* water = Particle::createWaterParticle(x * .1, y * .1 - 0.5 , z * .1);
        particle_grid->addToGrid(water);
      }
    }
  }

  for (int x = -60; x < 60; x++) {
    for (int z = -60; z < 60; z++) {
      Particle* boundary = Particle::createBoundaryParticle(x * .05, -1, z * .05);
      particle_grid->addToGrid(boundary);
    }
  }

  for (int y = -60; y < 60; y++) {
    for (int x = -60; x < 60; x++) {
      Particle* boundary = Particle::createBoundaryParticle(x * .05, y * 0.05, -3);
      particle_grid->addToGrid(boundary);
    }      
  }

  for (int y = -60; y < 60; y++) {
    for (int x = -60; x < 60; x++) {      
      Particle* boundary = Particle::createBoundaryParticle(x * 0.05, y * 0.05, 3);
      particle_grid->addToGrid(boundary);
    }
  }

  for (int z = -60; z < 60; z++) {
    for (int y = -60; y < 60; y++) {
      Particle* boundary = Particle::createBoundaryParticle(-3,  y * 0.05, z * 0.05);
      particle_grid->addToGrid(boundary);

    }
  }

  for (int z = -60; z < 60; z++) {
    for (int y = -60; y < 60; y++) {
      Particle* boundary = Particle::createBoundaryParticle(3, y * 0.05, z * 0.05);
      particle_grid->addToGrid(boundary);
    }      
  }

  //Calculate densities for particles
  #pragma omp parallel for
  for (int i = 0; i < particle_grid->water_particles->size(); ++i) {
    Particle* particle = particle_grid->water_particles->at(i);
    vector<Particle*>* neighbors = particle_grid->getNeighbors(particle);
    particle->set_density(neighbors);
    //delete neighbors;
  }

  //This initializes glut
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels. Add Depth Chanels
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

  // Initalize theviewport size
  viewport.w = 1920;
  viewport.h = 1200;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0,0);
  glutCreateWindow(argv[0]);

  initScene();              // quick function to set up scene

  glutReshapeFunc(myReshape);       // function to run when the window gets resized
  glutDisplayFunc(myDisplay);       // function to run when its time to draw something
  glutKeyboardFunc(myKeyboardFunc); // basic keys callback
  glutSpecialFunc(mySpecialKeyFunc); //special keys callback
  glutIdleFunc(myFrameMove);                   // function to run when not handling any other task

  glutMainLoop();             // infinite loop that will keep drawing and resizing
  // and whatever else

  return 0;
}








