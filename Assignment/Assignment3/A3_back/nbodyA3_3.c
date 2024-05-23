#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Define graviational constant and small number for smoother the simulation

#define G 6.67330e-11
//define as per given instruction
#define EPSILON 1e-3
#define DT 1e-5

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
	int N;
	//G = 100 / N;
	Vector dist_norm = distance_normalize(*p1, *p2);
	//Given G = 100/N
	double f = - (100/N) * ((p1->mass * p2->mass) / ((dist_fn + EPSILON) * (dist_fn + EPSILON) * (dist_fn + EPSILON)));
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
	//double dt = 0.01;	//Time step
	printf("\nPlease enter the Particle Number (N) and number of steps(nsteps): ");
	//printf("\nPlease enter the number of steps(nsteps): ");

	scanf("%d %d", &N, &nsteps);
	//scanf("%d",  &nsteps);
	
	//Read initial positions, velocities and mass of particles from file
	FILE *read_file = fopen("input_data/ellipse_N_00100.gal", "rb");
	if(read_file == NULL) {
		perror("Error in opening file!\n");
		return 1;
	}
	//fread(&N, sizeof(int), 1, read_file);
	Particle *particles = malloc(N * sizeof(Particle));	//Allocate Particle memory location
	if(particles == NULL) {
		perror("Memory allocation failed.\n");
		fclose(read_file);
		free(particles);
		return 1;
	}
	int num_pt = 0;
	for(int i = 0; i < N; i++) {
		fscanf(read_file, "%lf %lf %lf %lf %lf %lf", &particles[i].position.x, &particles[i].position.y, &particles[i].mass, &particles[i].velocity.x, &particles[i].velocity.y, &particles[i].brightness);
		//fread(&particles[i], sizeof(Particle), N, read_file);
		//fread(&particles[i], sizeof(Particle), 1, read_file);
		num_pt += 1;

	}
	
	
	//int particles_read = fread(particles, sizeof(Particle), N, read_file);
	fclose(read_file);
	
	if(num_pt != N) {
		fprintf(stderr, "Error: Read %d particles, expected %d particles.\n", num_pt, N);
		free(particles);
		
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
			particles[i].velocity.x += DT * particles[i].acceleration.x;
			particles[i].velocity.y += DT * particles[i].acceleration.y;
			
			//Update Potision
			particles[i].position.x += DT * particles[i].velocity.x;
			particles[i].position.y += DT * particles[i].velocity.y;
		}
		//Output of particle positions
		//for (int i = 0; i < N; i++) {
		//	printf("%f %f\n", particles[i].position.x, particles[i].position.y);
		//}
		//printf("\n");	
	}
	//Output in the result file
	FILE *output_file = fopen("ref_output_data/ellipse_result_N_00100.gal", "wb");
	if(output_file == NULL) {
		perror("Error in opening the file!\n");
		free(particles);
		fclose(output_file);
		return 1;
	}
	//double brightness = 0.0;

	for( int i = 0; i < N; i++) {
		//fprintf(output_file, "%lf %lf %lf %lf %lf\n", particles[i].mass, particles[i].position.x, particles[i].position.y, particles[i].velocity.x, particles[i].velocity.y);
		fprintf(output_file, "%lf %lf %lf %lf %lf %lf\n", particles[i].position.x, particles[i].position.y, particles[i].mass, particles[i].velocity.x, particles[i].velocity.y, particles[i].brightness);

		//fwrite(&particles[i], sizeof(Particle), 1, output_file);
	}
	
	//fwrite(particles, sizeof(Particle), N, output_file);
	fclose(output_file);
	
	free(particles);	//Free allocated memory
	return 0;
}

