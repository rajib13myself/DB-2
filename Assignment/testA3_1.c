#include <stdio.h>
#include <math.h>

#define N 100 // Maximum number of particles
#define G 6.67430e-11 // Gravitational constant
#define EPSILON 1e-3 // Small constant to avoid singularity

// Particle structure
typedef struct {
    double x;
    double y;
    double mass;
} Particle;

// Function to calculate distance between particles
double distance(Particle p1, Particle p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Function to calculate force exerted on particle i by other particles
void calculate_force(Particle particles[], int num_particles, double forces[][2]) {
    for (int i = 0; i < num_particles; i++) {
        double force_x = 0.0;
        double force_y = 0.0;
        for (int j = 0; j < num_particles; j++) {
            if (j != i) {
                double rij = distance(particles[i], particles[j]);
                double rij_normalized_x = (particles[i].x - particles[j].x) / rij;
                double rij_normalized_y = (particles[i].y - particles[j].y) / rij;
                double rij_epsilon = rij + EPSILON;
                double force_magnitude = -G * particles[i].mass * particles[j].mass / pow(rij_epsilon, 3);
                force_x += force_magnitude * rij_normalized_x;
                force_y += force_magnitude * rij_normalized_y;
            }
        }
        forces[i][0] = force_x;
        forces[i][1] = force_y;
    }
}

int main() {
    Particle particles[N];
    double forces[N][2];

    // Initialize particle positions and masses
    // (This is just an example, you should initialize these values according to your problem)
    for (int i = 0; i < N; i++) {
        particles[i].x = i;
        particles[i].y = i;
        particles[i].mass = 1.0;
    }

    // Calculate forces
    calculate_force(particles, N, forces);

    // Print forces
    for (int i = 0; i < N; i++) {
        printf("Force on particle %d: (%f, %f)\n", i, forces[i][0], forces[i][1]);
    }

    return 0;
}

