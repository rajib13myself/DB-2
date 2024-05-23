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

// Function to compute the gravitational force between two particles
void calculate_force(Particle p1, Particle p2, double *force_x, double *force_y) {
    double G = 6.67430e-11; // Gravitational constant
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double distance_squared = dx * dx + dy * dy;
    double distance = sqrt(distance_squared);
    double force_magnitude = G * p1.mass * p2.mass / (distance_squared);
    *force_x = force_magnitude * (dx / distance);
    *force_y = force_magnitude * (dy / distance);
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

    // Allocate memory for force matrix
    double **force_matrix = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        force_matrix[i] = (double*)malloc(2 * N * sizeof(double));
    }

    // Compute the force matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double force_x, force_y;
            if (i != j) {
                calculate_force(particles[i], particles[j], &force_x, &force_y);
                force_matrix[i][j] = force_x;
                force_matrix[i][N + j] = force_y;
            } else {
                force_matrix[i][j] = 0.0; // No force on itself
                force_matrix[i][N + j] = 0.0;
            }
        }
    }

    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        // Compute the acceleration using force matrix
        double *acceleration_x = (double*)calloc(N, sizeof(double));
        double *acceleration_y = (double*)calloc(N, sizeof(double));
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 2 * N; j++) {
                acceleration_x[i] += force_matrix[i][j] / particles[i].mass;
                acceleration_y[i] += force_matrix[i][N + j] / particles[i].mass;
            }
        }

        // Update velocities and positions
        for (int i = 0; i < N; i++) {
            particles[i].vx += acceleration_x[i] * delta_t;
            particles[i].vy += acceleration_y[i] * delta_t;
            particles[i].x += particles[i].vx * delta_t;
            particles[i].y += particles[i].vy * delta_t;
        }

        // Free allocated memory for acceleration
        free(acceleration_x);
        free(acceleration_y);
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
    for (int i = 0; i < N; i++) {
        free(force_matrix[i]);
    }
    free(force_matrix);
    free(particles);

    return 0;
}
