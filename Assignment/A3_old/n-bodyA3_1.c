#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Define gravity constant and small number for solution stability

#define G 6.67330e-11
#define EPSILON 1e-3

//Define particle position, velocity, acceleration and mass

typedef struct {
	double x, y;		//Particle position
	double velx, vely;	//Particle Velocity
	double accx, accy;	//Particle acceleration
	double mass;		//Particle mass
}Particle;

int N;			//Number of Particles

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
	p1->accx += F * (p2->x - p1->x);
	p2->accy += F * (p2->y - p2->y);
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

int main() {
	int N = 100;			//Define Number of particles
	
	Particle *particles = malloc(N*sizeof(Particle));		//Allocate memory location for particles.
	if(particles ==NULL) {
		perror("Memory allocation failed.\n");
		return 1;
	}	
	/*
	particles[0] = (Particle){ 0, 0, 0, 0, 0, 1000};
	particles[1] = (Particle){ 1, 0, 0, 0, 0, 100};
	*/
	//Initialize particles
	int k = 1000;	//Particles initialize points
	for (int i = 0; i < N; i++) {
		particles[i] = (Particle){ i, 0, 0, 0, 0, k};
		k -= 10;
	
	}
	
	//Simulation Parameters and accept user input value for entry N, dt and nsteps
	double dt = 0.01;		//Time Step
	int nsteps = 100;		//Number of steps
	/*
	printf("\nPlease Enter Two(2) values for time steps(dt) and number of steps(nsteps) with single space: ");
	scanf("%lf %d\n",  &dt, &nsteps);

	/*if(scanf("%lf %d\n",  &dt, &nsteps) != 2) {
		//Handle input error;
		printf("Error: Input ofrmat is incorrect!\n");
		return 1;
	}*/

	/*(if(dt > 1.00) {
		printf("Entered value for dt is higher value which is beyond this process.\n");
		return 1;
	}*/	

	
	
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

	//final position status
	for( int i = 0; i < N; i++) {
		printf("Particle %d: Position (%0.2f, %0.2f)\n", i, particles[i].x, particles[i].y);
	}

	//Free the allocated memory
	free(particles);
	return 0;
}


	
		
	


