#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>

#define G 100.0 // gravitational constant adjusted by number of particles
#define EPSILON 1e-3 // small number for stability

pthread_mutex_t mutex;

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
    double delta_t;
} ThreadArgs;

double distance(Particle *p1, Particle *p2) {
    return sqrt(((p2->x - p1->x) * (p2->x - p1->x))  + ((p2->y - p1->y) * (p2->y - p1->y)));
}

double force(Particle *p1, Particle *p2, int N) {
    double r = distance(p1, p2);
    double r_e = r + EPSILON;
    return -(100.0/N) * p1->mass * p2->mass /(r_e * r_e * r_e);
}

void* simulation_force(void* args) {
    ThreadArgs* threadArgs = (ThreadArgs*)args;
    Particle* particles = threadArgs->particles;
    int start = threadArgs->start;
    int end = threadArgs->end;
    double delta_t = threadArgs->delta_t; // Extract delta_t from ThreadArgs
    int N = threadArgs->N;
    double f_ij;
    double ax, ay;
    for (int i = start; i < end; i ++) {
        ax = 0.0; ay = 0.0;
        for (int j = 0; j < N; j++) {
            double rx = particles[i].x - particles[j].x;
            double ry = particles[i].y - particles[j].y;
            if ( j != i) {
                f_ij = force(&particles[i], &particles[j],N);
                ax += f_ij * rx/particles[i].mass;
                ay += f_ij * ry/particles[i].mass;
            }
        }
        particles[i].vx += ax * delta_t;
        particles[i].vy += ay * delta_t;
    }
    return NULL;
}

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

    pthread_mutex_init(&mutex, NULL);

    Particle *particles = (Particle*)malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 1;
    }

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

    #pragma omp parallel for num_threads(num_threads)
    for (int step = 0; step < nsteps; step++) {
        #pragma omp parallel for
        for(int i = 0; i < N; i++) {
            ThreadArgs args;
            args.particles = particles;
            args.N = N;
            args.start = 0;
            args.end = N;
            args.delta_t = delta_t; // Assign delta_t to ThreadArgs
            simulation_force((void*)&args);
        }
      }
      #pragma omp barrier
      #pragma omp parallel for
      for (int i = 0; i < N; i++) {
         particles[i].x += (particles[i].vx * delta_t);
         particles[i].y += (particles[i].vy * delta_t);
      }

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

    free(particles);
    pthread_mutex_destroy(&mutex);

    return 0;
}

