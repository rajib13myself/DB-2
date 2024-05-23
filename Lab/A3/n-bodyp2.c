#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100 // Number of particles
#define G (100.0 / N) // Gravitational constant scaled inversely with N
#define EPSILON 0.001 // Smoothing factor
#define DT 0.00001 // Time step size

typedef struct {
    double x, y; // Position
    double vx, vy; // Velocity
    double ax, ay; // Acceleration
    double mass; // Mass
    double brightness; //brightness
} Particle;

Particle particles[N];

void initialize_particles() {
    FILE *input_file;
    input_file = fopen("input_data/ellipse_N_00010.gal", "rb");
    if (input_file == NULL) {
        fprintf(stderr, "Error opening input file!\n");
        exit(1);
    }
    
    for (int i = 0; i < N; i++) {
        fscanf(input_file, "%lf %lf %lf %lf %lf %lf", &particles[i].x, &particles[i].y, &particles[i].mass,
                                                  &particles[i].vx, &particles[i].vy,
                                                  &particles[i].brightness);
    }
    
    fclose(input_file);
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

void write_results() {
    FILE *output_file;
    output_file = fopen("output.gal", "wb");
    if (output_file == NULL) {
        fprintf(stderr, "Error opening output file!\n");
        exit(1);
    }
    
    for (int i = 0; i < N; i++) {
        fprintf(output_file, "%lf %lf %lf %lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].mass,
                                                      particles[i].vx, particles[i].vy,
                                                      particles[i].brightness);
    }
    
    fclose(output_file);
}

int main() {
    initialize_particles();
    
    for (int step = 0; step < 10000; step++) {
        calculate_forces();
        update_particles();
    }
    
    write_results();
    
    return 0;
}

