#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


//Define graviational constant and small number for smoother the simulation

//#define G 6.67330e-11
#define G 100.0
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
Vector vector_cal(Particle p1, Particle p2) {
	Vector r1_p1p2;
	//r1_p1p2 = ((p1.x - p2.x) * uv.ex) + ((p1.y - p2.y) * uv.ey)
	//r1_p1p2.x = (p2.position.x - p1.position.x) * uv.ex;
	//r1_p1p2.y = (p2.position.y - p1.position.y) * uv.ey;
	r1_p1p2.x = (p2.position.x - p1.position.x);
	r1_p1p2.y = (p2.position.y - p1.position.y);

	return r1_p1p2;
}

//fucntion of distance between two particles
double distance_cal(Particle p1, Particle p2) {
	double dx = p2.position.x - p1.position.x;
	double dy = p2.position.y - p1.position.y;
	double dist_val= sqrt((dx * dx) + (dy * dy));
	if(dist_val <= 0) {
		return 0.0001;	//return small value
	}
	return (dist_val);
}

//function of normalize distance calculation
Vector distance_normalize(Vector r_ij) {
	//double dist_fn = distance_cal(p1, p2);
	double dist_fn = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y);
	Vector r2_p1p2 = {0};
	//UnitVector uv_xy = {1.0, 1.0};
	//Vector vect_fn = vector_cal(p1, p2);
	if (dist_fn != 0) {
		r2_p1p2.x = r_ij.x / dist_fn;
		r2_p1p2.y = r_ij.y / dist_fn;
	}
	return r2_p1p2;
}

//function to allocate memory for particle

Particle* allocate_particles(int N) {
	Particle* particles = malloc(N * sizeof(Particle));
	if(particles == NULL) {
		fprintf(stderr, "Memory allocation failed.\n");
		exit(1);
	}
	return particles;
}

//function to read particle from file 

void read_particles(const char* filename, Particle* particles, int N) {
	FILE* read_file = fopen(filename, "rb");
	if (read_file == NULL) {
		fprintf(stderr, "Error to opening the file.\n");
		exit(1);
	}
	for (int i = 0; i < N; i++) {
	
		
		//size_t  read_filedata = fread(&particles[i].position.x, sizeof(double), 1, read_file);
		
		fread(&particles[i].position.x, sizeof(double), 1, read_file);
		fread(&particles[i].position.y, sizeof(double), 1, read_file);
		fread(&particles[i].mass, sizeof(double), 1, read_file);
		fread(&particles[i].velocity.x, sizeof(double), 1, read_file);
		fread(&particles[i].velocity.y, sizeof(double), 1, read_file);
		fread(&particles[i].brightness, sizeof(double), 1, read_file);
		
		//read_filedata += fread(&particles[i].position.y, sizeof(double), 1, read_file);
		//read_filedata += fread(&particles[i].mass, sizeof(double), 1, read_file);
		//read_filedata += fread(&particles[i].velocity.x, sizeof(double), 1, read_file);
		//read_filedata += fread(&particles[i].velocity.y, sizeof(double), 1, read_file);
		//read_filedata += fread(&particles[i].brightness, sizeof(double), 1, read_file);
		
		//fread(&particles[i], sizeof(Particle), 1, read_file);
		//fscanf(read_file, "%lf %lf %lf %lf %lf %lf\n", &particles[i].position.x, &particles[i].position.y, &particles[i].mass, &particles[i].velocity.x, &particles[i].velocity.y, &particles[i].brightness);
	}
	
	fclose(read_file);
}



//function of force calculation of the particles
Vector force_particle(Particle *p1, Particle *p2) {
	double dist_fn = distance_cal(*p1, *p2);
	Vector r_ij = vector_cal(*p1, *p2);
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
	Vector dist_norm = distance_normalize(r_ij);
	//Given G = 100/N
	double f = - (G/N) * ((p1->mass * p2->mass) / ((dist_fn + EPSILON) * (dist_fn + EPSILON) * (dist_fn + EPSILON)));
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

int main(int argc, const char* args[]) {
	//int N, nsteps;		//Particle Number and number of steps variable
	//double dt = 0.01;	//Time step
	//printf("\nPlease enter the Particle Number (N) and number of steps(nsteps): ");
	if( argc != 6) {
		printf("Please enter argument parameters as 'N filename nsteps delta_t graphics(0/1)'\n");
		return 1;
	}
	int N = atoi(args[1]);
	const char* filename = args[2];
	int nsteps = atoi(args[3]);
	double delta_t = atof(args[4]);
	int graphics = atoi(args[5]);
	/*
	FILE *read_file = fopen(filename, "rb");
	if( read_file == NULL) {
		printf("Memory allocation failed \n");
		return 1;
	}

	Particle *particles = malloc(N * sizeof(Particle));	//Allocate Particle memory location
	if(particles == NULL) {
		perror("Memory allocation failed.\n");
		fclose(read_file);
		free(particles);
		return 1;
	}
	*/
	Particle* particles = allocate_particles(N);
	/*
	int particles_read = fread(particles, sizeof(Particle), N, read_file);
	fclose(read_file);
	if(particles_read != N) {
		printf("Error: read %d particles, expected %d particles.\n", particles_read, N);	
		free(particles);
		return 1;
	}
 	*/
	read_particles(filename, particles, N);

	//printf("\nPlease enter the number of steps(nsteps): ");

	//scanf("%d %d", &N, &nsteps);
	//scanf("%d",  &nsteps);
	
	//Read initial positions, velocities and mass of particles from file
	/*
	FILE *read_file = fopen("input_data/ellipse_N_00010.gal", "rb");
	if(read_file == NULL) {
		perror("Error in opening file!\n");
		return 1;
	}
	//fread(&N, sizeof(int), 1, read_file);
	*/
	/*
	int num_pt = 0;
	for(int i = 0; i < N; i++) {
		//fscanf(read_file, "%lf %lf %lf %lf %lf", &particles[i].mass, &particles[i].position.x, &particles[i].position.y, &particles[i].velocity.x, &particles[i].velocity.y);
		//fread(&particles[i], sizeof(Particle), N, read_file);
		fread(&particles[i], sizeof(Particle), 1, read_file);
		num_pt += 1;

	}
	
	
	
	//int particles_read = fread(particles, sizeof(Particle), N, read_file);
	fclose(read_file);
	
	if(num_pt != N) {
		fprintf(stderr, "Error: Read %d particles, expected %d particles.\n", num_pt, N);
		free(particles);
		
		return 1;
	}

	*/	
	
	
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
		//}
		//for (int i = 0; i < N; i++) {
			//Update Acceleration
			if(particles[i].mass !=0) {
				particles[i].acceleration.x = total_force.x / particles[i].mass;
				particles[i].acceleration.y = total_force.y / particles[i].mass;
			} else {
				particles[i].acceleration.x = 0;
				particles[i].acceleration.y = 0;
			}
						
			//Update Velocity
			particles[i].velocity.x += delta_t * particles[i].acceleration.x;
			particles[i].velocity.y += delta_t * particles[i].acceleration.y;
			
			//Update Potision
			particles[i].position.x += delta_t * particles[i].velocity.x;
			particles[i].position.y += delta_t * particles[i].velocity.y;
			//rounded the position values
			//particles[i].position.x = round(particles[i].position.x * 1e6) / 1e6;
			//particles[i].position.y = round(particles[i].position.y * 1e6) / 1e6;

		}
		//Output of particle positions
		//for (int i = 0; i < N; i++) {
		//	printf("%f %f\n", particles[i].position.x, particles[i].position.y);
		//}
		//printf("\n");	
	}
	//Output in the result file
	FILE *output_file = fopen("result.gal", "wb");
	if(output_file == NULL) {
		perror("Error in opening the file!\n");
		free(particles);
		fclose(output_file);
		return 1;
	}

	/*
	// Calculate the total size of data to write in bytes
	size_t data_size = N * sizeof(double) * 6;

    	// Allocate memory for the data
	double* data_buffer = malloc(data_size);
    	if (data_buffer == NULL) {
		perror("Memory allocation failed.\n");
	        free(particles);
	        fclose(output_file);
	        return 1;
	}

	// Copy the data to the buffer
	for (int i = 0; i < N; i++) {
		memcpy(data_buffer + i * 6, &particles[i].position.x, sizeof(double));
	        memcpy(data_buffer + i * 6 + 1, &particles[i].position.y, sizeof(double));
	        memcpy(data_buffer + i * 6 + 2, &particles[i].mass, sizeof(double));
	        memcpy(data_buffer + i * 6 + 3, &particles[i].velocity.x, sizeof(double));
	        memcpy(data_buffer + i * 6 + 4, &particles[i].velocity.y, sizeof(double));
	        memcpy(data_buffer + i * 6 + 5, &particles[i].brightness, sizeof(double));
	}

	// Write the data buffer to the file
    	fwrite(data_buffer, sizeof(double), N * 6, output_file);
	*/

    	// Free the allocated memory
    	//free(data_buffer);

	

	
	for( int i = 0; i < N; i++) {
		//fprintf(output_file, "%08lf %0lf %lf %08lf %08lf %lf\n",  particles[i].position.x, particles[i].position.y, particles[i].mass, particles[i].velocity.x, particles[i].velocity.y, particles[i].brightness);
		fwrite(&particles[i].position.x, sizeof(double), 1, output_file);
		fwrite(&particles[i].position.y, sizeof(double), 1, output_file);
		fwrite(&particles[i].mass, sizeof(double), 1, output_file);
		fwrite(&particles[i].velocity.x, sizeof(double), 1, output_file);
		fwrite(&particles[i].velocity.y, sizeof(double), 1, output_file);
		fwrite(&particles[i].brightness, sizeof(double), 1, output_file);

		//write each particle position(x,y), mass, velocity(x,y), brightness in this format
		//fprintf(output_file, "%08lf %08lf ", particles[i].position.x, particles[i].position.y);
		
	}
	
	//fwrite(particles, sizeof(Particle), N, output_file);
	fclose(output_file);
	
	free(particles);	//Free allocated memory
	return 0;
}

