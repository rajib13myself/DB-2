#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

        // Calculate forces between particles
        for (int i = 0; i < N; i++) {
            forces_x[i] = 0.0;
            forces_y[i] = 0.0;
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    double rx = particles[j].x - particles[i].x;
                    double ry = particles[j].y - particles[i].y;
                    double r = sqrt(rx * rx + ry * ry) + 1e-3;
                    double f_ij = (100.0 / N) * particles[i].mass * particles[j].mass / (r * r * r);
                    forces_x[i] += f_ij * rx / particles[i].mass;
                    forces_y[i] += f_ij * ry / particles[i].mass;
                }
            }
        }

        // Calculate acceleration and update velocities
        for (int i = 0; i < N; i++) {
            double ax = forces_x[i];
            double ay = forces_y[i];
            particles[i].vx += ax * delta_t;
            particles[i].vy += ay * delta_t;
        }

        // Free memory for forces
        free(forces_x);
        free(forces_y);

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
    free(particles);

    return 0;
}
