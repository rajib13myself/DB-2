#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <GL/glut.h>

//Define gravity constant and small number for solution stability

//#define G 6.67330e-11
#define G 100.0 / N		//Gravity constant scale by the number of particle as per assuignment Task-2
#define EPSILON 1e-3		//Small number for stability
#define L 1.0			//Domain Size as per given instruction
#define W 1.0			//

//Define particle position, velocity, acceleration and mass

typedef struct {
	double x, y;		//Particle position
	double velx, vely;	//Particle Velocity
	double accx, accy;	//Particle acceleration
	double mass;		//Particle mass
}Particle;

int N;			//Number of Particles
Particle *particles;

//Prototype of functions
double distance_func(Particle p1, Particle p2);
void fouce_particle(Particle *p1, Particle *p2);
void symplectic_method(Particle *p, double dt);
void drawScreen();
void updateScreen(int val);

//Define Screen dimention
const int WINDOW_WIDTH = 600;
const int WINDOW_HEIGHT = 400;



//define distance function for 2 particles
double distance_func(Particle p1, Particle p2) {
	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	double r_ij = sqrt(dx*dx + dy*dy);
	return r_ij;
}

//Calculate force component for a particle
void force_particle(Particle *p1, Particle *p2 ) {
	double r1_ij = distance_func(*p1, *p2);
	double F = -G * (p1->mass * p2->mass) / ((r1_ij + EPSILON)*(r1_ij + EPSILON)*(r1_ij + EPSILON));	//Stability solution with Plummer Spheres.
	p1->accx += F * (p2->x - p1->x);		//Particle P1 in x and y both directions
	p1->accy += F * (p2->y - p1->y);
	p2->accx -= F * (p2->x - p1->x);		//Particle P2 in both x and y directions
	p2->accy -= F * (p2->y - p2->y);
}

//Update velocity, positions by using Symplectic Euler method
void symplectic_method(Particle *p, double dt) {
	p->velx += dt * p->accx;		//Modify velocity by acceleration at x position
	p->vely += dt * p->accy;		//Modify velocity by acceleration at y position
	
	p->x += dt * p->velx;			//Modify x position by Velocity
	p->y += dt * p->vely;			//Modify y position by velocity

	p->accx = 0;				//Reset Accleration
	p->accy = 0;				//Reset Acceleration
}

//draw the graphics output screen
void drawScreen() {
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 1.0f, 1.0f); 		//Draw particle with white color
	glBegin(GL_POINTS);
	for (int i = 0; i < N; i++) {
		glVertex2f(particles[i].x, particles[i].y);
	}
	glEnd(); 
	
	glFlush();
}

//Graphics screen update function
void updateScreen(int val) {
	//update particles position and velocity
	for(int ts = 0; ts < 1; ts++) {
		for (int i = 0; i < N; i++) {
			for( int j = 0; j < N; j++) {
				if(j != i) {
					force_particle(&particles[i], &particles[j]);
				}
			}
		}
		for(int i =0; i < N; i++) {
			symplectic_method(&particles[i], 1e-5);
		}
	}
	glutPostRedisplay();
	glutTimerFunc(25, updateScreen, 0);
}

	
int main(int argc, char *argv[]) {
	if( argc !=6) {
		printf("Uses: %s N filename nsteps delta_t graphics\n", argv[0]);
		return 1;
	}

	//int N = 100;			//Define Number of particles
	//printf("Enter value for Number of particles(N): ");
	//scanf("%d\n", &N);
	//five input initialization
	N = atoi(argv[1]);				//Number of particles
	char *filename = malloc(strlen(argv[2]) +1);	//Filename from initial conditions from
	strcpy(filename, argv[2]);
	int nsteps = atoi(argv[3]);			//Number of timesteps
	double dt = atof(argv[4]);			//Timestep size
	int graphics = atoi(argv[5]);			//Flag for graphics output 
	
	if (N <= 0) {
		printf("Error: Number of particle (N) must be greater than 0\n");
		return 1;
	}
	
	Particle *particles = malloc(N* sizeof(Particle));		//Allocate memory location for particles.
	if(particles ==NULL) {
		perror("Memory allocation failed.\n");
		return 1;
	}	
	/*
	particles[0] = (Particle){ 0, 0, 0, 0, 0, 1000};
	particles[1] = (Particle){ 1, 0, 0, 0, 0, 100};
	*/
	//Initialize particles from file
	FILE *input_file = fopen(filename, "rb");
	if (input_file == NULL) {
		perror("Error: input file unable to find in the location.\n");
		free(particles);
		free(filename);
		return 1;
	}
	for(int i = 0; i < N; i++) {
		//fscanf(input_file, "%lf %lf %lf %lf %lf", &particles[i].x, &particles[i].y, &particles[i].velx, &particles[i].vely, &particles[i].mass);
		fread(&particles[i], sizeof(Particle), 1, input_file);
	}
	fclose(input_file);
	
		
	
	//Simulation Parameters and accept user input value for entry N, dt and nsteps
	//double dt = 1e-5;		//Time Step
	//int nsteps = 100;		//Number of steps
	/*
	printf("\nPlease Enter Two(2) values for time steps(dt) and number of steps(nsteps) with single space: ");
	scanf("%lf %d\n",  &dt, &nsteps);

	/*if(scanf("%lf %d\n",  &dt, &nsteps) != 2) {
		//Handle input error;
		printf("Error: Input ofrmat is incorrect!\n");
		return 1;
	}*/

	if(dt > 1.00) {
		printf("Entered value for dt is higher value which is beyond this process.\n");
		return 1;
	}	

	
	
	//Simulation process
	for(int ts = 0; ts < nsteps; ts++) {
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				if(j != i) {
					force_particle(&particles[i], &particles[j]);
				}
			}
		}
	
	//Modify velocity and position by Using the symplectic Euler time integration method
		for( int i = 0; i < N; i++) {
		symplectic_method(&particles[i], dt);
		}
	}
	
	
	//Innitialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutCreateWindow("Galaxy Simulation");
	
	//Setting OpenGL projection Matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, L, 0.0, W, -1.0, 1.0);
	
	//Setting OpenGL Modelview Matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	//Screen Background Color
	glClearColor(0.0f, 0.0f,0.0f, 1.0f);

	//Setting OpenGL callback functions
	glutDisplayFunc(drawScreen);
	glutTimerFunc(25, updateScreen, 0);

	//GLUT main loop start
	glutMainLoop();
	
	//Output in file with give instruction.
	FILE *output_file = fopen("result.gal", "w");
	if(output_file == NULL) {
		perror("Error opening output in the result file.\n");
		free(particles);
		free(filename);
		return 1;
	}
	
	for(int i = 0; i < N; i++) {
		//fprintf(output_file, "Particle %d: Position (%0.2f, %0.2f), Velocity (%0.2f, %0.2f), Mass %0.2f\n", i, particles[i].x, particles[i].y, particles[i].velx, particles[i].vely, particles[i].mass);
		fwrite(&particles[i], sizeof(Particle), 1, output_file);
	}
	
	/*//final position status
	for( int i = 0; i < N; i++) {
		printf("Particle %d: Position (%0.2f, %0.2f)\n", i, particles[i].x, particles[i].y);
	}*/
	//CLose file
	fclose(output_file);
	
	//Free the allocated memory
	free(particles);
	free(filename);
	return 0;
}


	
		
	


