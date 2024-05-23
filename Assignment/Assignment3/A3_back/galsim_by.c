#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define G 100.0 // gravitational constant adjusted by number of particles
//#define EPSILON 1e-3 // small number for stability


typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    double brightness; // brightness
} Particle;

// Function to read particles from binary file
int read_particles(const char *filename, Particle *particles, int num_particles) {
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        return 0;
    }

    int success = fread(particles, sizeof(Particle), num_particles, input_file);
    fclose(input_file);

    return success;
}

// Function to write particles to binary file
int write_particles(const char *filename, Particle *particles, int num_particles) {
    FILE *output_file = fopen(filename, "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        return 0;
    }

    int success = fwrite(particles, sizeof(Particle), num_particles, output_file);
    fclose(output_file);

    return success;
}



int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
        return 0;
    }
    
    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    
    // Allocate memory for particles
    Particle *particles = (Particle*)malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 0;
    }
    /*
    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        free(particles);
        return 0;
    }
    	
    for (int i = 0; i < N ; i++) {
        fread(&particles[i].x, sizeof(double), 1, input_file);
        fread(&particles[i].y, sizeof(double), 1, input_file);
        fread(&particles[i].mass, sizeof(double), 1, input_file);
        fread(&particles[i].vx, sizeof(double), 1, input_file);
        fread(&particles[i].vy, sizeof(double), 1, input_file);
        fread(&particles[i].brightness, sizeof(double), 1, input_file);
        
        // Print values for debugging
        //printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, particles[i].brightness);
        //fprintf(input_file, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf", &particles[i].x, &particles[i].y; &particles[i].mass, &particles[i].vx, &particles[i].vy, &particles[i].brightness);
        
    }*/
    // Read initial conditions from file
    int success_read = read_particles(filename, particles, N);
    if (!success_read) {
        free(particles);
        return 0;
    }
        
    //fclose(input_file);
    
    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        double f_ij;
        double ax, ay;
        for (int i = 0; i < N; i ++) {
            //double f_ij;
	    ax = 0.0;  ay = 0.0;
            for (int j = 0; j < N; j++) {
                double rx=particles[i].x-particles[j].x;
                double ry=particles[i].y-particles[j].y;
		double r = sqrt(rx * rx + ry * ry);
                if ( j != i) {
                    f_ij = -(100.0/N) * (particles[i].mass) * particles[j].mass /((r + 1e-3) *  (r + 1e-3) * (r + 1e-3));
                    ax += f_ij * rx/particles[i].mass;
                    ay += f_ij * ry/particles[i].mass;
                }
                
            }
            particles[i].vx += ax * delta_t;
            particles[i].vy += ay * delta_t;
           // printf("%f   %f\n",ax,ay);
           // double dummy;
            //scanf("stop:",&dummy);
        }
        for (int i = 0; i < N; i++) {
            // Update positions
            particles[i].x += (particles[i].vx * delta_t);
            particles[i].y += (particles[i].vy * delta_t);
            
        }
    }
    /*
    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
	fclose(output_file);
        free(particles);
        return 0;
    }
    
    for (int i = 0; i < N; i++) {
        fwrite(&particles[i].x, sizeof(double), 1, output_file);
        fwrite(&particles[i].y, sizeof(double), 1, output_file);
        fwrite(&particles[i].mass, sizeof(double), 1, output_file);
        fwrite(&particles[i].vx, sizeof(double), 1, output_file);
        fwrite(&particles[i].vy, sizeof(double), 1, output_file);
        fwrite(&particles[i].brightness, sizeof(double), 1, output_file);
        
        // Print values for debugging
        //printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, particles[i].brightness);
        
    }*/
    // Write results to file
    int success_write = write_particles("result.gal", particles, N);
    if (!success_write) {
        printf("Error: Unable to write output file.\n");
        free(particles);
        return 0;
    }
        
    //fclose(output_file);
    
    // Free allocated memory
    free(particles);
    
    return 0;
}

