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
    
    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
	fclose(input_file);
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
        
    }
        
    fclose(input_file);
    
    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        register double f_ij;
        register double ax, ay;
	register double r, rx, ry;
        for (int i = 0; i < N; i ++) {
            ax = 0.0;  ay = 0.0;
            for (int j = 0; j < N; j++) {
                rx = (particles[j].x-particles[i].x);
		ry = (particles[j].y-particles[i].y);
		//Debug value of x, y
                printf("rx and ry for i,j: %lf , %lf\n", rx, ry);
		r = sqrt(rx * rx + ry * ry) + 1e-3;
                if ( j != i) {
                    f_ij = (100.0/N) * (particles[i].mass) * particles[j].mass /(r * r * r);
                    ax += f_ij * rx/particles[i].mass;
                    ay += f_ij * ry/particles[i].mass;
                }
                
            }
            particles[i].vx += ax * delta_t;
            particles[i].vy += ay * delta_t;
            printf("%f   %f\n",ax,ay);
           // double dummy;
            //scanf("stop:",&dummy);
        }
        for (int i = 0; i < N; i++) {
            // Update positions
            particles[i].x += (particles[i].vx * delta_t);
            particles[i].y += (particles[i].vy * delta_t);
       }
    }
    
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
        
    }
       
    fclose(output_file);
    
    // Free allocated memory
    free(particles);
    
    return 0;
}

