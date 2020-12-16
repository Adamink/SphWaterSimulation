#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <omp.h>

#include "constants.h"
#include "sph.h"

static SPH::Sph* sph_instance;

static void initGLutWindow(){
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);

	using namespace visual_const;
	gluPerspective(fovy, aspect, zFar, zNear);
	gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz);
}

static void display(){
	initGLutWindow();
	sph_instance->step();
	if(sph_instance->getStep() > program_const::kTotalStep){
		exit(0);
	}
	glutSwapBuffers();
}

static void idle(){
	glutPostRedisplay();
}

static void runWithNoVisualization(){
	while(sph_instance->getStep() <= program_const::kTotalStep){
		sph_instance->step();
	}
}
int main(int argc, char* argv[]){
	using namespace program_const;
	using namespace visual_const;
	omp_set_num_threads(NUM_THREADS_COMPUTING);
	sph_instance = new SPH::Sph();
	if(IF_VISUALIZE){
		glutInit(&argc, argv);
		glutInitWindowSize(kInitWindowSizeX, kInitWindowSizeY);
		glutInitWindowPosition(kInitWindowPositionX, kInitWindowPositionY);
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);
		glutCreateWindow("SphWaterSimulation");
		initGLutWindow();
		glutDisplayFunc(display);
		glutIdleFunc(idle);
		glutMainLoop();
	}
	else{
		runWithNoVisualization();
	}
	return 0;
}

