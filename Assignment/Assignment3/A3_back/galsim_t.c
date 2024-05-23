#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 100.0 // gravitational constant adjusted by number of particles
#define EPSILON 1e-3 // small number for stability


typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    //    double ax, ay; //acceleration
    double brightness; // brightness
} Particle;

// Function to calculate distance between two particles
double distance(Particle *p1, Particle *p2) {
    return sqrt(((p2->x - p1->x) * (p2->x - p1->x))  + ((p2->y - p1->y) * (p2->y - p1->y)));
    //return sqrt(((p2.x - p1.x) * (p2.x - p1.x))  + ((p2.y - p1.y) * (p2.y - p1.y)));
    //return sqrt((pow(p2->x - p1->x), 2)  + (pow(p2->y - p1->y), 2));
    // return sqrt(((p1.x - p2.x) * (p1.x - p2.x))  + ((p1.y - p2.y) * (p1.y - p2.y)));
}

//Function of distance vector
double dist_vec(Particle *p1, Particle *p2) {
    //return ((p2->x - p1->x) * 1.0 + (p2->y - p1->y)) * 1.0;		//ex & ey as unique vector
    return ((p2->x - p1->x)  + (p2->y - p1->y)) ;		//ex & ey as unique vector
    //return ((p1.x - p2.x) * 1.0 + (p1.y - p2.y)) * 1.0;		//ex & ey as unique vector
}


// Function to calculate total force acting on a particle
double force(Particle *p1, Particle *p2,int N) {
    double r = distance(p1, p2);
    return -(100.0/N) * p1->mass * p2->mass /pow(r + EPSILON, 3);
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
        
        // Print values for debugging
        printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, particles[i].brightness);
        //fprintf(input_file, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf", &particles[i].x, &particles[i].y; &particles[i].mass, &particles[i].vx, &particles[i].vy, &particles[i].brightness);
        
    }
    
    fclose(input_file);
    
    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        //calculate_force(particles, N, delta_t);
        double f_ij;
        //F_ij = 0;
        double ax = 0, ay = 0;
        for (int i = 0; i < N; i ++) {
            ax = 0.0; ay = 0.0;
            for (int j = 0; j < N; j++) {
                double rx=particles[i].x-particles[j].x;
                double ry=particles[i].y-particles[j].y;
                if ( j != i) {
                    f_ij = force(&particles[i], &particles[j],N);
                    ax += f_ij * rx/particles[i].mass;
                    ay += f_ij * ry/particles[i].mass;
                }
                
            }
            particles[i].vx += ax * delta_t;
            particles[i].vy += ay * delta_t;
            printf("%f   %f\n",ax,ay);
            double dummy;
            //scanf("stop:",&dummy);
        }
        for (int i = 0; i < N; i++) {
            // Update positions
            particles[i].x += (particles[i].vx * delta_t);
            particles[i].y += (particles[i].vy * delta_t);
            
        }
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
        
        // Print values for debugging
        printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, particles[i].brightness);
        
    }
    
    fclose(output_file);
    
    // Free allocated memory
    free(particles);
    
    return 0;
}
