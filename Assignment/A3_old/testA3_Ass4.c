#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.67430e-11

typedef struct {
    double x;
    double y;
    double mass;
    double vx;
    double vy;
    double brightness;
} Particle;

void update_positions(Particle *particles, int N, double dt) {
    for (int i = 0; i < N; i++) {
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
    }
}

void update_velocities(Particle *particles, int N, double dt) {
    for (int i = 0; i < N; i++) {
        double ax = 0, ay = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double r = sqrt(dx*dx + dy*dy);
                double F = G * particles[i].mass * particles[j].mass / (r*r);
                ax += F * dx / r;
                ay += F * dy / r;
            }
        }
        particles[i].vx += ax * dt;
        particles[i].vy += ay * dt;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage: ./galsim N filename nsteps delta_t graphics\n");
        return 1;
    }
    
    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double dt = atof(argv[4]);
    int graphics = atoi(argv[5]);
    
    Particle *particles = malloc(N * sizeof(Particle));
    if (!particles) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    // Read initial configuration from file
    FILE *file = fopen(filename, "rb");
    if (!file) {
        printf("Failed to open file %s\n", filename);
        free(particles);
        return 1;
    }
    fread(particles, sizeof(Particle), N, file);
    fclose(file);
    
    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        update_positions(particles, N, dt);
        update_velocities(particles, N, dt);
    }
    
    // Save final configuration to result.gal
    file = fopen("result.gal", "wb");
    if (!file) {
        printf("Failed to open result file for writing\n");
        free(particles);
        return 1;
    }
    fwrite(particles, sizeof(Particle), N, file);
    fclose(file);
    
    free(particles);
    
    return 0;
}

