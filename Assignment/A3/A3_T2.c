#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 100.0 // gravitational constant adjusted by number of particles
#define EPSILON 1e-3 // small number for stability
#define DT 1e-5 // time step size

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double ax, ay; // acceleration
    double mass; // mass
    double brightness; //brightness
} Particle;


// Function to calculate distance between two particles
double distance(Particle p1, Particle p2) {
    return sqrt(((p1.x - p2.x) * (p1.x - p2.x))  + ((p1.y - p2.y) * (p1.y - p2.y)));
}

//Function for distance vector(Particle p1, Particle p2)
double dist_vector(Particle p1, Particle p2) {
    return ((p1.x - p2.x) * 1.0 + (p1.y - p2.y) * 1.0); //Here ex and ey are unit vector
}

//Function for normalize both particles distance
double norm_dist(Particle p1, Particle p2) {
    double r_vec = dist_vector(p1, p2);
    double r_dist = distance(p1, p2);
    double r_hat = r_vec / r_dist;
    return r_hat;
}
// Function to calculate force between two particles using modified force formula
double force(Particle p1, Particle p2) {
    double r_dist = distance(p1, p2);
    double r_vec = dist_vector(p1, p2);
    double f = -G * p1.mass * p2.mass * ( r_vec / ((r_dist + EPSILON) * (r_dist + EPSILON) * (r_dist + EPSILON)));
    return f;
}


int main() {
    int N = 10;
    double L = 1.0, W = 1.0;
    
    /*
    // Read input parameters from file
    FILE *input_file = fopen("input_data/ellipse_N_00010.gal", "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open input file.\n");
        return 1;
    }
    fscanf(input_file, "%lf %lf", &N, &L, &W);
    fclose(input_file);
    */
    // Allocate memory for particles
    Particle *particles = (Particle*)malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 1;
    }
    
    // Read initial conditions from file
    FILE *input_file = fopen("input_data/ellipse_N_00010.gal", "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        free(particles);
        return 1;
    }
    for (int i = 0; i < N; i++) {
        fscanf(input_file, "%lf %lf %lf %lf %lf %lf", &particles[i].x, &particles[i].y, &particles[i].vx, &particles[i].vy, &particles[i].mass, &particles[i].brightness);
    }
    fclose(input_file);
    
    // Time loop
    for (double t = 0; t < 1.0; t += DT) {
        // Calculate forces and update accelerations
        for (int i = 0; i < N; i++) {
            double f;
            particles[i].ax = 0;
            particles[i].ay = 0;
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    double f = f + force(particles[i], particles[j]);
                    //double r = distance(particles[i], particles[j]);
                   // particles[i].ax += f * (particles[j].x - particles[i].x) / r;
                    //particles[i].ay += f * (particles[j].y - particles[i].y) / r;
                }
            }
            particles[i].ax += f * particles[i].mass;
            particles[i].ay += f * particles[i].mass; 
        }
        
        // Update velocities and positions
        for (int i = 0; i < N; i++) {
            if (i < N-1) {
                particles[i+1].vx = particles[i].vx + particles[i].ax * DT;    // / particles[i].mass;
                particles[i+1].vy = particles[i].vy + particles[i].ay * DT;    // particles[i].mass;
                particles[i+1].x = particles[i].x + particles[i+1].vx * DT;
                particles[i+1].y = particles[i].y + particles[i+1].vy * DT;
            }
        }
    }
    
    // Write results to file
    FILE *output_file = fopen("results.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        free(particles);
        return 1;
    }
    for (int i = 0; i < N; i++) {
        fprintf(output_file, "%lf %lf %lf %lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, particles[i].brightness);
    }
    fclose(output_file);
    
    // Free allocated memory
    free(particles);
    
    return 0;
}
