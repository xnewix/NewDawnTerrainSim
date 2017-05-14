#include <time.h>
#include "SpaceCamera.h"
#include "Grid.h"
#include "hgtreader.h"

bool mousein=false;
bool polygonmode=false;

CTexture texture_grass;

class Settings {
private:
	std::string grid_filename;
	double sunx, suny, sunz;
	double camerax, cameray, cameraz;
	bool lighting;
	
public:
	Settings() {
			setsettings();
	}
	~Settings() {}

	std::string getGrid_Filename() { 
		return grid_filename;
	}
	double getSunX() { return sunx; }
	double getSunY() { return suny; }
	double getSunZ() { return sunz; }

	double getCameraX() { return camerax; }
	double getCameraY() { return cameray; }
	double getCameraZ() { return cameraz; }

	bool   getLighting() { return lighting; }

	void setsettings() {
		std::ifstream datafile;
		std::string tmp;
		datafile.open("settings.txt");
		datafile >> grid_filename;
		datafile >> sunx;
		datafile >> suny;
		datafile >> sunz;

		datafile >> tmp;
		if (tmp.compare("lightingON")) { 
			lighting = true;
			tmp.clear();
		}
		else {
			lighting = false;
			tmp.clear();
		}

		datafile.close();
		tmp.clear();
	}
};

Settings worldsettings;

SpaceCamera spacecamera(30.0,90.0,30.0,-15.0,-135.0);
//30,90,30

CGrid cgrid1(513); 

int gridmap1;

//MapOSM maposm("map.osm");

void init() {
	glClearColor(0.0,0.0,0.0,0.5);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45,640.0/480.0,0.1,360000.0);  
	glMatrixMode(GL_MODELVIEW);

	glEnable(GL_DEPTH_TEST);	
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);	
	// glEnable(GL_TEXTURE_2D);

    GLfloat mat_specular[]   = {0.1, 0.1, 0.1, 1.0};
    GLfloat mat_diffuse[]    = {0.5, 0.5, 0.5, 1.0};
	GLfloat mat_ambient[]    = {0.5, 0.5, 0.5, 1.0};
	GLfloat mat_shininess[]  = {0.0};	

	GLfloat light_ambient[]  = {0.3, 0.3, 0.3, 1.0};
	//GLfloat light_ambient[]  = {0.0, 0.0, 0.0, 1.0};

	GLfloat light_diffuse[]  = {0.3, 0.3, 0.3, 1.0};
	//GLfloat light_diffuse[]  = {0.0, 0.0, 0.0, 1.0};

	GLfloat light_specular[] = {0.1, 0.1, 0.1, 1.0};
	//GLfloat light_specular[] = {0.0, 0.0, 0.0, 1.0};

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	glShadeModel(GL_SMOOTH);	

	texture_grass.createtexture("grass.bmp");
	texture_grass.usetexture();

	if (worldsettings.getLighting() == true) {
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
	}
    

	gridmap1 = cgrid1.draw(spacecamera);	
		

}

void display(SDL_Window *window)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	spacecamera.Control(window,0.2,0.2,mousein);	
	spacecamera.UpdateCamera();	

	if(polygonmode) {
		    glPolygonMode(GL_FRONT, GL_LINE);
            glPolygonMode(GL_BACK, GL_LINE);
	}
	if(!polygonmode) {
			glPolygonMode(GL_FRONT, GL_FILL);
            glPolygonMode(GL_BACK, GL_FILL);
	}

	GLfloat light_position[] = { worldsettings.getSunX(), worldsettings.getSunY(), worldsettings.getSunZ(), 1.0 }; 
	// light position was previously 0, 64, 0

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glCullFace(GL_BACK);
	glTranslated(0.0,0.0,0.0);
	glPushMatrix();
		glCallList(gridmap1);	
		//maposm.draw();	
	glPopMatrix();
}


int main(int argc, char **argv) {

    SDL_Init(SDL_INIT_VIDEO);   

	
	SDL_Window *window = SDL_CreateWindow("SDL2/OpenGL Demo", 
		SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 
		640, 480, 
		SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN 
	);  
	
	// | SDL_WINDOW_FULLSCREEN_DESKTOP

	SDL_GLContext glcontext = SDL_GL_CreateContext(window);

	bool running=true;
	Uint32 start;
	SDL_Event event;
	init();	
 
	OpenglVersion oglversion;
	oglversion.openglversion_output();
	
	while(running){
		start=SDL_GetTicks();
		while(SDL_PollEvent(&event)) {
			switch(event.type) {				
				case SDL_QUIT:
					running=false;
					break;
				case SDL_MOUSEBUTTONDOWN:
					mousein=true;
					SDL_ShowCursor(SDL_DISABLE);
					break;
				case SDL_KEYDOWN:					
					if(event.key.keysym.sym==SDLK_p) {
						mousein=false;
						SDL_ShowCursor(SDL_ENABLE);
						break;
					}
					if(event.key.keysym.sym==SDLK_ESCAPE) {
						running=false;
						break;
					}
					if((event.key.keysym.sym==SDLK_m) && (polygonmode==false)) {
						polygonmode = true;
						break;
					}
					if((event.key.keysym.sym==SDLK_m) && (polygonmode==true)) {
						polygonmode = false;
						break;
					}
			}
		}
		display(window);		
	    SDL_GL_SwapWindow(window);  
	}
	SDL_Quit();	
	return 0;	
}
