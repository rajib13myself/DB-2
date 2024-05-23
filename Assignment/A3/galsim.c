#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 100.0 // gravitational constant adjusted by number of particles
#define EPSILON 1e-3 // small number for stability

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    double brightness; // brightness
} Particle;

// Function to calculate distance between two particles
double distance(Particle p1, Particle p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}
double ax, ay;
// Function to calculate force between two particles using modified force formula
void calculate_force(Particle *particles, int N, double delta_t) {
    for (int i = 0; i < N; i++) {
        //ax = 0;
        //ay = 0;
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double r = distance(particles[i], particles[j]);
                double f = -G * particles[i].mass * particles[j].mass / pow(r + EPSILON, 3);
                ax += f * (particles[j].x - particles[i].x);
                ay += f * (particles[j].y - particles[i].y);
            }
        }
        particles[i].vx += ax * delta_t;
        particles[i].vy += ay * delta_t;
        particles[i].x += particles[i].vx * delta_t;
        particles[i].y += particles[i].vy * delta_t;
    }
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
    
    Particle *particles = (Particle*)malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 1;
    }
    
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open input file.\n");
        free(particles);
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        fread(&particles[i], sizeof(Particle), 1, input_file);
    }
    fclose(input_file);
    
    for (int step = 0; step < nsteps; step++) {
        calculate_force(particles, N, delta_t);
        
        if (graphics) {
            // Code to display graphics
        }
    }
    
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        free(particles);
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        fwrite(&particles[i], sizeof(Particle), 1, output_file);
    }
    fclose(output_file);
    
    free(particles);
    
    return 0;
}
