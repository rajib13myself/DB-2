#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100 // Number of particles
#define G 9.81 // Gravitational constant
#define EPSILON 0.001 // Smoothing factor
#define DT 0.01 // Time step size

typedef struct {
    double x, y; // Position
    double vx, vy; // Velocity
    double ax, ay; // Acceleration
    double mass; // Mass
} Particle;

Particle particles[N];

void initialize_particles() {
    for (int i = 0; i < N; i++) {
        particles[i].x = rand() % 100; // Random initial position
        particles[i].y = rand() % 100;
        particles[i].vx = 0; // Initial velocity is zero
        particles[i].vy = 0;
        particles[i].mass = 1; // Mass of each particle is 1 for simplicity
    }
}

void calculate_forces() {
    for (int i = 0; i < N; i++) {
        particles[i].ax = 0;
        particles[i].ay = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                // Calculate distance between particles
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double distance_squared = dx * dx + dy * dy;
                double distance = sqrt(distance_squared);
                
                // Calculate force
                double force_magnitude = (G * particles[i].mass * particles[j].mass) / (distance_squared + EPSILON);
                double fx = force_magnitude * dx / distance;
                double fy = force_magnitude * dy / distance;
                
                // Update acceleration
                particles[i].ax += fx / particles[i].mass;
                particles[i].ay += fy / particles[i].mass;
            }
        }
    }
}

void update_particles() {
    for (int i = 0; i < N; i++) {
        // Update velocity
        particles[i].vx += particles[i].ax * DT;
        particles[i].vy += particles[i].ay * DT;
        
        // Update position
        particles[i].x += particles[i].vx * DT;
        particles[i].y += particles[i].vy * DT;
    }
}

int main() {
    initialize_particles();
    
    for (int step = 0; step < 1000; step++) {
        calculate_forces();
        update_particles();
        
        // Output particle positions (optional)
        if (step % 100 == 0) {
            printf("Step %d:\n", step);
            for (int i = 0; i < N; i++) {
                printf("Particle %d: x=%f, y=%f\n", i, particles[i].x, particles[i].y);
            }
            printf("\n");
        }
    }
    
    return 0;
}

