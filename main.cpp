
#include <vector>
#include <queue>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

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

//CONSTANTS
long double PI = atan(1)*4;
long double E = 2.7182818284590452353;
long double WATER_MASS = 1.0;
long double WATER_VICOSITY_COEFFICIENT = 1.0;
long double WATER_BUOYANCY_STRENGTH = 1.0;
long double WATER_GAS_CONSTANT = 1.0;
long double WATER_REST_DENSITY = 1.0;
long double FOG_MASS = 1.0;
long double FOG_VICOSITY_COEFFICIENT = 1.0;
long double FOG_BUOYANCY_STRENGTH = 1.0;
long double FOG_GAS_CONSTANT = 1.0;
long double FOG_REST_DENSITY = 1.0;

using namespace std;

//The vertices that define our polygon
vector<vector<pair<ThreeDVector*, ThreeDVector*> > > polygons;

//Save Boolean
bool save = false;
//Filename
static const char* file_name;

//Max/Min x,y,z
long double max_x = numeric_limits<long double>::min();
long double min_x = numeric_limits<long double>::max();
long double max_y = numeric_limits<long double>::min();
long double min_y = numeric_limits<long double>::max();
long double max_z = numeric_limits<long double>::min();
long double min_z = numeric_limits<long double>::max();

//Print Function for debugging
void print(string _string) {
  cout << _string << endl;
}

void print(long double num) {
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
        long double x = atof(splitline[1].c_str());
        long double y = atof(splitline[2].c_str());
        long double z = atof(splitline[3].c_str());
        ThreeDVector* vertex = new ThreeDVector(x, y, z);
        vertices.push_back(vertex);
      }
      //vn x y z
      else if(!splitline[0].compare("vn")) {
        long double x = atof(splitline[1].c_str());
        long double y = atof(splitline[2].c_str());
        long double z = atof(splitline[3].c_str());
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
          char *pch = std::strtok(v_str, "/");
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

  long double radius = max(max_x - min_x, max(max_y - min_y, max_z - min_x)) / 2;
  long double center_x = (max_x + min_x) / 2;
  long double center_y = (max_y + min_y) / 2;
  long double center_z = (max_z + min_z) / 2;

  long double multiplier = 1;

  //glOrtho((center_x - radius) * multiplier, (center_x + radius) * multiplier, (center_y - radius) * multiplier, (center_y + radius) * multiplier, 1.0, 1.0 + radius);
  gluPerspective(90, float(w)/float(h), 1.0, 1.0 + radius * 10);
  glViewport (0,0,viewport.w,viewport.h);

}

//****************************************************
// Simple init function
//****************************************************
void initScene(){
  glClearColor(0.52941f, 0.80784f, 0.98039f, 1.0f); // Clear to black, fully transparent

  // Enable lighting and the light we have set up
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glEnable(GL_LIGHT1);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_RESCALE_NORMAL);

  //Set lighting parameters
  GLfloat light_position0[] = {0.0, 1.0, -1.0, 0};
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

  //Set Material Parameters
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

  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  myReshape(viewport.w,viewport.h);
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);       // clear the color buffer

  glMatrixMode(GL_MODELVIEW);         // indicate we are specifying camera transformations
  glLoadIdentity();                   // make sure transformation is "zero'd"

  long double radius = max(max_x - min_x, max(max_y - min_y, max_z - min_x)) / 2;
  long double center_x = (max_x + min_x) / 2;
  long double center_y = (max_y + min_y) / 2;
  long double center_z = (max_z + min_z) / 2;

  //gluLookAt(center_x + 1000000, center_y, center_z - 2500000, center_x + 200000, center_y, center_z, 0, 1, 0);

  gluLookAt(220, 0, -150, 175, 0, 0, 0, 1, 0);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Enable shading
  glShadeModel(GL_SMOOTH);

  long double scale_factor = 1.0/10000.0;
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

  if (save) {
    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    vector<unsigned char> buf(w * h * 4);

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, &buf[0] );

    createPng(file_name, buf, w, h);
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

int main(int argc, char *argv[]) {

  for (int i = 10; i < argc; i++) {
    if (string(argv[i]) == "-save") {
      if(i + 1 < argc) {
        save = true;
        file_name = argv[i + 1];
        i = i + 1;
      }
    } 
  }

  //Parse Polygons the Golden Gate
  parseObj("Golden Gate Bridge.obj");

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

  glutMainLoop();             // infinite loop that will keep drawing and resizing
  // and whatever else

  return 0;
}








