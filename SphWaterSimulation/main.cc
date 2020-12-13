#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <omp.h>

#include "constants.h"
#include "sph.h"

SPH::Sph* sph0 = new SPH::Sph();
#define NUM_THREADS 4

bool Compute_code = true;
bool Read_file = !Compute_code;

void init(void){
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);

	using namespace camera_const;
	gluPerspective(fovy, aspect, zFar, zNear);
	gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz);
}

static void display(void){
	init();

	if(Compute_code){ sph0->step(); }
	else if(Read_file){ sph0->OnlyReadFileStep(); }

	if(sph0->getstep() > sph0->getMaxstep()){
		system("pause");
		exit(0);
	}

	glutSwapBuffers();
}

static void idle(void){
	glutPostRedisplay();
}

int main(int argc, char* argv[]){
	omp_set_num_threads(NUM_THREADS);
	glutInit(&argc, argv);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(10, 10);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);

	glutCreateWindow("SPH method");

	init();

	glutDisplayFunc(display);
	glutIdleFunc(idle);

	glutMainLoop();

	return EXIT_SUCCESS;
}

