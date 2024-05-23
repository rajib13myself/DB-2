#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALIGNMENT 32 // Adjust alignment for aligned memory allocation

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    double brightness; // brightness
} Particle;

// Matrix-vector multiplication for computing forces
void matvec_force(double **mat_a, double *vec_b, double *vec_c, int SIZE) {
    for (int i = 0; i < SIZE; i++) {
        vec_c[i] = 0.0; // Reset force accumulator
        for (int j = 0; j < SIZE; j++) {
            vec_c[i] += mat_a[i][j] * vec_b[j];
        }
    }
}

// Matrix initialization
void init_matvec(double **mat_a, double *vec_b, double *vec_c, double *vec_d, int SIZE, Particle *particles) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (i != j) {
                double rx = particles[j].x - particles[i].x;
                double ry = particles[j].y - particles[i].y;
                double r = sqrt(rx * rx + ry * ry) + 1e-3;
                mat_a[i][j] = (100.0 / SIZE) * (particles[i].mass) * particles[j].mass / (r * r * r);
            } else {
                mat_a[i][j] = 0.0; // Diagonal elements are 0
            }
        }
        vec_b[i] = particles[i].mass;
        vec_c[i] = 0.0;
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

    // Matrix initialization
    double **force_matrix = (double**)malloc(N * sizeof(double*));
    double *masses = (double*)malloc(N * sizeof(double));
    double *forces_x = (double*)malloc(N * sizeof(double));
    double *forces_y = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        force_matrix[i] = (double*)malloc(N * sizeof(double));
    }

    // Initialize force matrix and vectors
    init_matvec(force_matrix, masses, forces_x, forces_y, N, particles);

    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        // Compute forces using matrix-vector multiplication
        matvec_force(force_matrix, masses, forces_x, N);
        matvec_force(force_matrix, masses, forces_y, N);

        // Update velocities and positions
        for (int i = 0; i < N; i++) {
            particles[i].vx += forces_x[i] * delta_t;
            particles[i].vy += forces_y[i] * delta_t;
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
        for (int i = 0; i < N; i++) {
            free(force_matrix[i]);
        }
        free(force_matrix);
        free(masses);
        free(forces_x);
        free(forces_y);
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
    for (int i = 0; i < N; i++) {
        free(force_matrix[i]);
    }
    free(force_matrix);
    free(masses);
    free(forces_x);
    free(forces_y);

    return 0;
}
