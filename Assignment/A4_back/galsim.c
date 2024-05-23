#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>

#define G 100.0 // gravitational constant adjusted by number of particles
//#define EPSILON 1e-3 // small number for stability

//Define mutex for thread synchronization
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    double brightness; // brightness
} Particle;

typedef struct {
    int start;
    int end;
    int N;
    Particle *particles;
} ThreadArgs;

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Usage: %s N filename nsteps delta_t graphics n_thread\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    int num_threads = atoi(argv[6]);

    // Initialize mutex
    pthread_mutex_init(&mutex, NULL);

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
        fclose(input_file);
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

    // Parallelize Simulation loop using OpenMP
    #pragma omp parallel for num_threads(num_threads)
    for (int step = 0; step < nsteps; step++) {
        for(int i = 0; i < N; i++) {
            double ax = 0.0, ay = 0.0;
            for (int j = 0; j < N; j++) {
                double rx = particles[j].x - particles[i].x;
                double ry = particles[j].y - particles[i].y;
                double r = sqrt(rx * rx + ry * ry) + 1e-3;
                if ( j != i) {
                    double f_ij = (100.0/N) * particles[i].mass * particles[j].mass /(r * r * r);
                    ax += f_ij * rx / particles[i].mass;
                    ay += f_ij * ry / particles[i].mass;
                }
            }
            particles[i].vx += ax * delta_t;
            particles[i].vy += ay * delta_t;
        }

        // Update particles position
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < N; i++) {
            particles[i].x += (particles[i].vx * delta_t);
            particles[i].y += (particles[i].vy * delta_t);
        }
    }

    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        fclose(output_file);
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
    //Destroy mutex
    pthread_mutex_destroy(&mutex);

    return 0;
}

