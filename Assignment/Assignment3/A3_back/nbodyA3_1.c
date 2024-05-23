#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Define graviational constant and small number for smoother the simulation

#define G 6.67330e-11
#define EPSILON 1e-3

//define structure for 2D vector

typedef struct {
	double x, y;			//particle position
}Vector;

//define structure of a particle

typedef struct {
	double mass;
	Vector position;
	Vector velocity;
	Vector acceleration;
	double brightness;
}Particle;

typedef struct {
	float ex, ey;
}UnitVector;

//function of vector calculation
Vector vector_cal(Particle p1, Particle p2, UnitVector uv) {
	Vector r1_p1p2;
	//r1_p1p2 = ((p1.x - p2.x) * uv.ex) + ((p1.y - p2.y) * uv.ey)
	r1_p1p2.x = (p1.position.x - p2.position.x) * uv.ex;
	r1_p1p2.y = (p1.position.y - p2.position.y) * uv.ey;
	return r1_p1p2;
}

//fucntion of distance between two particles
double distance_cal(Particle p1, Particle p2) {
	double dx = p1.position.x - p2.position.x;
	double dy = p1.position.y - p2.position.y;
	double dist_val= (dx * dx) + (dy * dy);
	if(dist_val <= 0) {
		return 0.0001;	//return small value
	}
	return sqrt(dist_val);
}

//function of normalize distance calculation
Vector distance_normalize(Particle p1, Particle p2) {
	double dist_fn = distance_cal(p1, p2);
	Vector r2_p1p2 = {0};
	UnitVector uv_xy = {1.0, 1.0};
	Vector vect_fn = vector_cal(p1, p2, uv_xy);
	if (dist_fn != 0) {
		r2_p1p2.x = vect_fn.x / dist_fn;
		r2_p1p2.y = vect_fn.y / dist_fn;
	}
	return r2_p1p2;
}

//function of force calculation of the particles
Vector force_particle(Particle *p1, Particle *p2) {
	double dist_fn = distance_cal(*p1, *p2);
	if(dist_fn == 0) {
		Vector def_force = {0,0};
		return def_force;
	}
	if (p1->mass == 0 || p2->mass == 0) {
		Vector def_force = {0, 0};
		return def_force;
	}

	Vector dist_norm = distance_normalize(*p1, *p2);
	double f = -G * (p1->mass * p2->mass) / ((dist_fn + EPSILON) * (dist_fn + EPSILON) * (dist_fn + EPSILON));
	Vector f_p1p2 = {0};
	f_p1p2.x = f * dist_norm.x;
	f_p1p2.y = f * dist_norm.y;
	
	if(isnan(f_p1p2.x) || isnan(f_p1p2.y)) {
		Vector def_force = {0, 0};
		return def_force;
	}
	return f_p1p2;
}

//Simulation total proces by using Symplectic Euler updates.

int main() {
	int N, nsteps;		//Particle Number and number of steps variable
	double dt = 0.01;	//Time step
	printf("\nPlease enter the Particle Number (N) and number of steps(nsteps): ");
	scanf("%d %d", &N, &nsteps);
	Particle *particles = malloc(N * sizeof(Particle));	//Allocate Particle memory location
	if(particles == NULL) {
		perror("Memory allocation failed.\n");
		return 1;
	}
	//Simulation loop
	for(int ts = 0; ts < nsteps; ts++) {
		//Caculate forces on each particle
		for(int i = 0; i < N; i++) {
			Vector total_force = {0.0, 0.0};
			for (int j = 0; j < N; j++) {
				if(j != i) {
					Vector f_ij = force_particle(&particles[i], &particles[j]);
					total_force.x += f_ij.x;
					total_force.y += f_ij.y;
				}
			}
			//Update Acceleration
			if(particles[i].mass !=0) {
				particles[i].acceleration.x = total_force.x / particles[i].mass;
				particles[i].acceleration.y = total_force.y / particles[i].mass;
			} else {
				particles[i].acceleration.x = 0;
				particles[i].acceleration.y = 0;
			}
						
			//Update Velocity
			particles[i].velocity.x += dt * particles[i].acceleration.x;
			particles[i].velocity.y += dt * particles[i].acceleration.y;
			
			//Update Potision
			particles[i].position.x += dt * particles[i].velocity.x;
			particles[i].position.y += dt * particles[i].velocity.y;
		}
		//Output of particle positions
		for (int i = 0; i < N; i++) {
			printf("%f %f\n", particles[i].position.x, particles[i].position.y);
		}
		printf("\n");	
	}

	free(particles);	//Free allocated memory
	return 0;
}

