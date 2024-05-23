#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#define NUM_THREADS 4


// Structure to hold thread-specific data
typedef struct {
    int id;
    int N;
    double *x;
    double *y;
    double *mass;
    double *vx;
    double *vy;
    double *brightness;
    int nsteps;
    double delta_t;
} ThreadData;

// Global variables
pthread_barrier_t barrier;

// Function to calculate forces
void *calculate_forces(void *thread_args) {
    ThreadData *data = (ThreadData *)thread_args;

    int start = (data->id * data->N) / NUM_THREADS;
    int end = ((data->id + 1) * data->N) / NUM_THREADS;

    for (int step = 0; step < data->nsteps; step++) {
        double rx, ry, r, f_ij;
        double *forces_x = calloc(data->N, sizeof(double));
        double *forces_y = calloc(data->N, sizeof(double));

        if (forces_x == NULL || forces_y == NULL) {
            printf("Error: Memory allocation failed.\n");
            free(forces_x);
            free(forces_y);
            return NULL;
        }

        // Calculate forces
        for (int i = start; i < end; i++) {
            f_ij = (100.0 / data->N) * data->mass[i];
            for (int j = 0; j < data->N; j++) {
                if (j != i) {
                    rx = data->x[i] - data->x[j];
                    ry = data->y[i] - data->y[j];
                    r = sqrt(rx * rx + ry * ry) + 1e-3; // Add epsilon to avoid division by zero
                    double f_ji = - f_ij * data->mass[j] / (r * r * r);
                    forces_x[i] += f_ji * rx / data->mass[i];
                    forces_y[i] += f_ji * ry / data->mass[i];
                }
            }
        }

        // Update velocities
        for (int i = start; i < end; i++) {
            data->vx[i] += forces_x[i] * data->delta_t;
            data->vy[i] += forces_y[i] * data->delta_t;
        }

        // Wait for all threads to finish calculations
        pthread_barrier_wait(&barrier);

        // Update positions
        for (int i = start; i < end; i++) {
            data->x[i] += data->vx[i] * data->delta_t;
            data->y[i] += data->vy[i] * data->delta_t;
        }

        // Wait for all threads to finish updating positions
        pthread_barrier_wait(&barrier);
    }

    return NULL;
}

// Main function
int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Usage: %s N filename nsteps delta_t graphics nthreads\n", argv[0]);
        return 0;
    }

    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    int nthreads = atoi(argv[6]);

    // Allocate memory for particle properties
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *mass = malloc(N * sizeof(double));
    double *vx = malloc(N * sizeof(double));
    double *vy = malloc(N * sizeof(double));
    double *brightness = malloc(N * sizeof(double));

    if (x == NULL || y == NULL || mass == NULL || vx == NULL || vy == NULL || brightness == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 0;
    }

    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        free(x);
        free(y);
        free(mass);
        free(vx);
        free(vy);
        free(brightness);
        fclose(input_file);
        return 0;
    }

    // Read particle data
    for (int i = 0; i < N; i++) {
        fread(&x[i], sizeof(double), 1, input_file);
        fread(&y[i], sizeof(double), 1, input_file);
        fread(&mass[i], sizeof(double), 1, input_file);
        fread(&vx[i], sizeof(double), 1, input_file);
        fread(&vy[i], sizeof(double), 1, input_file);
        fread(&brightness[i], sizeof(double), 1, input_file);
    }

    fclose(input_file);

    // Initialize thread data
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];

    pthread_barrier_init(&barrier, NULL, NUM_THREADS);

    // Create threads
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].id = i;
        thread_data[i].N = N;
        thread_data[i].x = x;
        thread_data[i].y = y;
        thread_data[i].mass = mass;
        thread_data[i].vx = vx;
        thread_data[i].vy = vy;
        thread_data[i].brightness = brightness;
        thread_data[i].nsteps = nsteps;
        thread_data[i].delta_t = delta_t;
        pthread_create(&threads[i], NULL, calculate_forces, (void *)&thread_data[i]);
    }

    // Wait for all threads to finish
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        free(x);
        free(y);
        free(mass);
        free(vx);
        free(vy);
        free(brightness);
        fclose(output_file);
        return 0;
    }

    for (int i = 0; i < N; i++) {
        fwrite(&x[i], sizeof(double), 1, output_file);
        fwrite(&y[i], sizeof(double), 1, output_file);
        fwrite(&mass[i], sizeof(double), 1, output_file);
        fwrite(&vx[i], sizeof(double), 1, output_file);
        fwrite(&vy[i], sizeof(double), 1, output_file);
        fwrite(&brightness[i], sizeof(double), 1, output_file);
    }

    fclose(output_file);

    // Free allocated memory
    free(x);
    free(y);
    free(mass);
    free(vx);
    free(vy);
    free(brightness);

    return 0;
}

