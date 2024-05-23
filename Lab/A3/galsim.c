#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics.h"

const float circleRadius=0.025, circleColor=0;
const int windowWidth=800;

void keep_within_box(float* xA, float* yA) {
  if(*xA > 1)
    *xA = 0;
  if(*yA > 1)
    *yA = 0;
}

typedef struct {
    double x, y; // Position
    double vx, vy; // Velocity
    double ax, ay; //accleration
    double mass; // Mass
    double brightness; // Brightness
} Particle;

void initialize_particles(Particle particles[], char *filename, int N) {
    FILE *read_file;
    read_file = fopen(filename, "rb");
    if (read_file == NULL) {
        fprintf(stderr, "Error opening input file!\n");
        exit(1);
    }

    //fread(particles, sizeof(Particle), N, input_file);
	for (int i = 0; i < N; i++) {
	 fread(&particles[i].x, sizeof(double), 1, read_file);
         fread(&particles[i].y, sizeof(double), 1, read_file);
         fread(&particles[i].mass, sizeof(double), 1, read_file);
         fread(&particles[i].vx, sizeof(double), 1, read_file);
         fread(&particles[i].vy, sizeof(double), 1, read_file);
         fread(&particles[i].brightness, sizeof(double), 1, read_file);
	}
    fclose(read_file);
}

void calculate_forces(Particle particles[], int N, double G, double EPSILON) {
    for (int i = 0; i < N; i++) {
        particles[i].ax = 0;
        particles[i].ay = 0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double distance_squared = dx * dx + dy * dy;
                double distance = sqrt(distance_squared);
                
                double force_magnitude = (G * particles[i].mass * particles[j].mass) / (distance_squared + EPSILON);
                double fx = force_magnitude * dx / distance;
                double fy = force_magnitude * dy / distance;
                
                particles[i].ax += fx / particles[i].mass;
                particles[i].ay += fy / particles[i].mass;
            }
        }
    }
}

void update_particles(Particle particles[], int N, double DT) {
    for (int i = 0; i < N; i++) {
        particles[i].vx += particles[i].ax * DT;
        particles[i].vy += particles[i].ay * DT;
        
        particles[i].x += particles[i].vx * DT;
        particles[i].y += particles[i].vy * DT;
    }
}

void write_results(Particle particles[], char *filename, int N) {
    FILE *output_file;
    output_file = fopen(filename, "wb");
    if (output_file == NULL) {
        fprintf(stderr, "Error opening output file!\n");
        exit(1);
    }
    
   // fwrite(particles, sizeof(Particle), N, output_file);
    for (int i = 0; i < N; i++) {

    fwrite(&particles[i].x, sizeof(double), 1, output_file);
    fwrite(&particles[i].y, sizeof(double), 1, output_file);
    fwrite(&particles[i].mass, sizeof(double), 1, output_file);
    fwrite(&particles[i].vx, sizeof(double), 1, output_file);
    fwrite(&particles[i].vy, sizeof(double), 1, output_file);
    fwrite(&particles[i].brightness, sizeof(double), 1, output_file);
   }
    
    fclose(output_file);
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
        return 1;
    }
    
    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    
    Particle *particles = malloc(N * sizeof(Particle));
    
    initialize_particles(particles, filename, N);
    
    double G = 100.0 / N;
    double EPSILON = 0.001;
    
    for (int step = 0; step < nsteps; step++) {
        calculate_forces(particles, N, G, EPSILON);
        update_particles(particles, N, delta_t);
        
        if (graphics == 1) {
            // Code for visualization (not provided)
        }
    }
    
    write_results(particles, "result.gal", N);
    
    free(particles);
    
    return 0;
}

