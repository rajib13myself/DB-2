#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALIGNMENT 16 // Adjust alignment for aligned memory allocation
#define VECTOR_SIZE 16 // Adjust vector size for SIMD optimization

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    double brightness; // brightness
} Particle;

void update_forces_and_velocities(Particle *particles, double *forces_x, double *forces_y, int N, double delta_t) {
    // Compute forces
    for (int i = 0; i < N; i++) {
        forces_x[i] = 0.0;
        forces_y[i] = 0.0;

        // Calculate forces for particles in vectorized manner
        for (int j = 0; j < N; j += VECTOR_SIZE) {
            double rx[VECTOR_SIZE], ry[VECTOR_SIZE], r[VECTOR_SIZE], f_ij[VECTOR_SIZE];

            for (int k = 0; k < VECTOR_SIZE; k++) {
                int idx = j + k;
                if (idx < N && i != idx) {
                    rx[k] = particles[idx].x - particles[i].x;
                    ry[k] = particles[idx].y - particles[i].y;
                    r[k] = sqrt(rx[k] * rx[k] + ry[k] * ry[k]) + 1e-3;
                    f_ij[k] = (100.0 / N) * particles[i].mass * particles[idx].mass / (r[k] * r[k] * r[k]);
                    forces_x[i] += f_ij[k] * rx[k];
                    forces_y[i] += f_ij[k] * ry[k];
                }
            }
        }
    }

    // Update velocities
    for (int i = 0; i < N; i++) {
        particles[i].vx += forces_x[i] / particles[i].mass * delta_t;
        particles[i].vy += forces_y[i] / particles[i].mass * delta_t;
    }
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
    Particle *particles = (Particle*)aligned_alloc(ALIGNMENT, N * sizeof(Particle));
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

    for (int i = 0; i < N; i++) {
        fread(&particles[i].x, sizeof(double), 1, input_file);
        fread(&particles[i].y, sizeof(double), 1, input_file);
        fread(&particles[i].mass, sizeof(double), 1, input_file);
        fread(&particles[i].vx, sizeof(double), 1, input_file);
        fread(&particles[i].vy, sizeof(double), 1, input_file);
        fread(&particles[i].brightness, sizeof(double), 1, input_file);
    }

    fclose(input_file);

    // Allocate memory for forces
    double *forces_x = (double*)malloc(N * sizeof(double));
    double *forces_y = (double*)malloc(N * sizeof(double));
    if (forces_x == NULL || forces_y == NULL) {
        printf("Error: Memory allocation failed.\n");
        free(forces_x);
        free(forces_y);
        free(particles);
        return 0;
    }

    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        // Compute forces and update velocities
        update_forces_and_velocities(particles, forces_x, forces_y, N, delta_t);

        // Update positions
        for (int i = 0; i < N; i++) {
            particles[i].x += particles[i].vx * delta_t;
            particles[i].y += particles[i].vy * delta_t;
        }
    }

    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        fclose(output_file);
        free(forces_x);
        free(forces_y);
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
    }

    fclose(output_file);

    // Free allocated memory
    free(forces_x);
    free(forces_y);
    free(particles);

    return 0;
}
