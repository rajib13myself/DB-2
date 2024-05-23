#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define gravity constant and small number for solution stability
#define G 100.0 / N
#define EPSILON 1e-3
#define L 1.0
#define W 1.0

// Define particle structure
typedef struct {
    double x, y;
    double velx, vely;
    double accx, accy;
    double mass;
} Particle;

int main() {
  //Open input file to read number of particles
    char filename[] = "/home/ubuntu/HPP/Assignment/A3/circles_N_2.gal"; 	
    FILE *input_file = fopen(filename, "r");
    if (input_file == NULL) {
        perror("Error: opening input file");
        return 1;
    }

    int N;
    if (fscanf(input_file, "%d", &N) != 1) {
        printf("Error: reading number of particles from input file");
        fclose(input_file);
        return 1;
    }
    fclose(input_file);

    if (N <= 0) {
        printf("Error: Number of particles (N) must be greater than 0.\n");
        return 1;
    }

    // Allocate memory for particles
    Particle *particles = malloc(N*sizeof(Particle));
    if (particles == NULL) {
        perror("Memory allocation failed\n");
        return 1;
    }

    // Initialize particles
    double mass_step = 1000.0 / N;
    for (int i = 0; i < N; i++) 	{
        particles[i] = (Particle){i * (double)N * L, 0, 0, 0, 0, 0,(mass_step * (N - i))};
    }

    // Simulation parameters as per instruction
    double dt = 1e-5;
    int nsteps = 100;

    // Simulation loop
    for (int ts = 0; ts < nsteps; ts++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    double dx = particles[j].x - particles[i].x;
                    double dy = particles[j].y - particles[i].y;
                    double r_ij = sqrt(dx * dx + dy * dy);
                    double F = -G * (particles[i].mass * particles[j].mass)* ((r_ij + EPSILON) * (r_ij + EPSILON) * (r_ij + EPSILON));
                    particles[i].accx += F * dx;
                    particles[i].accy += F * dy;
                    particles[j].accx -= F * dx;
                    particles[j].accy -= F * dy;
                }
            }
        }

        for (int i = 0; i < N; i++) {
            particles[i].velx += dt * particles[i].accx;
            particles[i].vely += dt * particles[i].accy;
            particles[i].x += dt * particles[i].velx;
            particles[i].y += dt * particles[i].vely;
            particles[i].accx = 0;
            particles[i].accy = 0;
        }
    }

    // Write particle data to output file
    FILE *output_file = fopen("/home/Ubuntu/HPP/Assignment/A3/results.txt", "w");
    if (output_file == NULL) {
        perror("Error: opening output file.\n");
        free(particles);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fprintf(output_file, "Particle %d Position (%0.2f, %0.2f), Velocity (%0.2f, %0.2f), Mass %0.2f\n",
                i, particles[i].x, particles[i].y, particles[i].velx, particles[i].vely, particles[i].mass);
    }

    fclose(output_file);
    free(particles);

    return 0;
}
