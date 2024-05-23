#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 100.0 // gravitational constant adjusted by number of particles
#define EPSILON 1e-3 // small number for stability
#define TIMESTEP 0.1    // time step size

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    double brightness; // brightness
} Particle;

// Function to calculate distance between two particles
double distance(Particle *p1, Particle *p2) {
    return sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2));
}

// Function to calculate normalized distance vector
void normalized_distance_vector(Particle *p1, Particle *p2, double *nx, double *ny) {
    double r = distance(p1, p2);
    *nx = (p1->x - p2->x) / r;
    *ny = (p1->y - p2->y) / r;
}

// Function to calculate total force acting on a particle
double force(Particle *p1, Particle *p2) {
    double r = distance(p1, p2);
    int N;
    return -(G/N) * p1->mass * p2->mass * pow(r + EPSILON, -3);
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
        return 1;
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
        return 1;
    }

    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        free(particles);
        return 1;
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
        // Update forces and accelerations
	double f = 0.0;
        for (int i = 0; i < N; i++) {
            double ax = 0.0;
            double ay = 0.0;
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    f = force(&particles[i], &particles[j]);
                    double nx, ny;
                    normalized_distance_vector(&particles[i], &particles[j], &nx, &ny);
                    ax += f * nx / particles[i].mass;
                    ay += f * ny / particles[i].mass;
                }
            }
            // Update velocities
            particles[i].vx += ax * delta_t;
            particles[i].vy += ay * delta_t;

            // Update positions
            particles[i].x += particles[i].vx * delta_t;
            particles[i].y += particles[i].vy * delta_t;
        }
	
        // Output particle positions at each time step (optional)
        // for (int i = 0; i < N; i++) {
        //     printf("Particle %d: x = %f, y = %f\n", i, particles[i].x, particles[i].y);
        // }
    }

    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        free(particles);
        return 1;
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

