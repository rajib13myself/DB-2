#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double data[6]; // 1D array for position (x, y), velocity (vx, vy), mass, brightness
} Particle;

// Function to calculate forces with vectorization
void calculate_forces_vectorized(int N, Particle* particles, double* forces_x, double* forces_y) {
    for (int i = 0; i < N; i++) {
        forces_x[i] = 0.0;
        forces_y[i] = 0.0;

        for (int j = 0; j < i +1; j++) {
            if (j != i) {
                double rx = particles[j].data[0] - particles[i].data[0];
                double ry = particles[j].data[1] - particles[i].data[1];
                double r_squared = rx * rx + ry * ry + 1e-3; // Add epsilon to avoid division by zero
                double f_ij = (100.0 / N) * particles[i].data[4] * particles[j].data[4] / (r_squared * sqrt(r_squared));

                forces_x[i] += f_ij * rx;
                forces_y[i] += f_ij * ry;
            }
        }
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
    Particle *particles = (Particle*)malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 0;
    }

    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        free(particles);
        return 0;
    }

    for (int i = 0; i < N ; i+=6) {
        fread(&particles[i], sizeof(double), 1, input_file);
        fread(&particles[i+1], sizeof(double), 1, input_file);
        fread(&particles[i+2], sizeof(double), 1, input_file);
        fread(&particles[i+3], sizeof(double), 1, input_file);
        fread(&particles[i+4], sizeof(double), 1, input_file);
        fread(&particles[i+5], sizeof(double), 1, input_file);
        // Print values for debugging
        printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i], particles[i+1], particles[i+2], particles[i+3], particles[i+4], particles[i+5]);
    }

    fclose(input_file);

    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        // Allocate memory for forces
        double *forces_x = (double*)malloc(N * sizeof(double));
        double *forces_y = (double*)malloc(N * sizeof(double));
        if (forces_x == NULL || forces_y == NULL) {
            printf("Error: Memory allocation failed.\n");
            free(forces_x);
            free(forces_y);
            return 0;
        }

        // Calculate forces between particles using vectorization
        calculate_forces_vectorized(N, particles, forces_x, forces_y);

        // Calculate acceleration and update velocities
        for (int i = 0; i < N; i++) {
            particles[i].data[2] += forces_x[i] * delta_t;
            particles[i].data[3] += forces_y[i] * delta_t;
        }

        // Free memory for forces
        free(forces_x);
        free(forces_y);

        // Update positions
        for (int i = 0; i < N; i++) {
            particles[i].data[0] += particles[i].data[2] * delta_t;
            particles[i].data[1] += particles[i].data[3] * delta_t;
        }
    }

    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        free(particles);
        return 0;
    }

    for (int i = 0; i < N; i+=6) {
        fwrite(&particles[i], sizeof(double), 1, output_file);
        fwrite(&particles[i+1], sizeof(double), 1, output_file);
        fwrite(&particles[i+2], sizeof(double), 1, output_file);
        fwrite(&particles[i+3], sizeof(double), 1, output_file);
        fwrite(&particles[i+4], sizeof(double), 1, output_file);
        fwrite(&particles[i+5], sizeof(double), 1, output_file);
        
        // Print values for debugging
        printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i], particles[i+1], particles[i+2], particles[i+3], particles[i+4], particles[i+5]);
    }

    fclose(output_file);

    // Free allocated memory
    free(particles);

    return 0;
}
