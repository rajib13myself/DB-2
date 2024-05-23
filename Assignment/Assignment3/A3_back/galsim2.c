#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Define gravitational constant and small number for smoothing the simulation
#define G 6.67330e-11
#define EPSILON 1e-3
#define DT 1e-5

// Define structure for 2D vector
typedef struct {
    double x, y; // Particle position
} Vector;

// Define structure of a particle
typedef struct {
    double mass;
    Vector position;
    Vector velocity;
    Vector acceleration;
    double brightness;
} Particle;

// Define structure for unit vector
typedef struct {
    float ex, ey;
} UnitVector;

// Function for vector calculation
Vector vector_cal(Particle p1, Particle p2, UnitVector uv) {
    Vector r1_p1p2;
    r1_p1p2.x = (p2.position.x - p1.position.x) * uv.ex;
    r1_p1p2.y = (p2.position.y - p1.position.y) * uv.ey;
    return r1_p1p2;
}

// Function for calculating distance between two particles
double distance_cal(Particle p1, Particle p2) {
    double dx = p2.position.x - p1.position.x;
    double dy = p2.position.y - p1.position.y;
    double dist_val = (dx * dx) + (dy * dy);
    return sqrt(dist_val);
}

// Function for normalizing distance calculation
Vector distance_normalize(Particle p1, Particle p2) {
    double dist_fn = distance_cal(p1, p2);
    Vector r2_p1p2 = {0};
    UnitVector uv_xy = {1.0, 1.0};
    Vector vect_fn = vector_cal(p1, p2, uv_xy);
    if (dist_fn != 0) {
        r2_p1p2.x = vect_fn.x / dist_fn;
        r2_p1p2.y = vect_fn.y / dist_fn;
    }
    return r2_p1p2;
}

// Function to allocate memory for particles
Particle* allocate_particles(int N) {
    Particle* particles = malloc(N * sizeof(Particle));
    if (particles == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    return particles;
}

// Function to read particles from file
void read_particles(const char* filename, Particle* particles, int N) {
    FILE* read_file = fopen(filename, "rb");
    if (read_file == NULL) {
        fprintf(stderr, "Error opening the file.\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        size_t read_filedata = fread(&particles[i].position.x, sizeof(double), 1, read_file);
        read_filedata += fread(&particles[i].position.y, sizeof(double), 1, read_file);
        read_filedata += fread(&particles[i].mass, sizeof(double), 1, read_file);
        read_filedata += fread(&particles[i].velocity.x, sizeof(double), 1, read_file);
        read_filedata += fread(&particles[i].velocity.y, sizeof(double), 1, read_file);
        read_filedata += fread(&particles[i].brightness, sizeof(double), 1, read_file);
    }
    fclose(read_file);
}

// Function for calculating force between particles
Vector force_particle(Particle *p1, Particle *p2, int N) {
    double dist_fn = distance_cal(*p1, *p2);
    if (dist_fn == 0 || p1->mass == 0 || p2->mass == 0) {
        Vector def_force = {0, 0};
        return def_force;
    }
    Vector dist_norm = distance_normalize(*p1, *p2);
    double f = - (100.0 / N) * ((p1->mass * p2->mass) / pow((dist_fn + EPSILON), 3));
    Vector f_p1p2 = {0};
    f_p1p2.x = f * dist_norm.x;
    f_p1p2.y = f * dist_norm.y;
    return f_p1p2;
}

int main(int argc, const char* args[]) {
    if (argc != 6) {
        printf("Please enter argument parameters as 'N filename nsteps delta_t graphics(0/1)'\n");
        return 1;
    }
    int N = atoi(args[1]);
    const char* filename = args[2];
    int nsteps = atoi(args[3]);
    double delta_t = atof(args[4]);
    int graphics = atoi(args[5]);

    Particle* particles = allocate_particles(N);
    read_particles(filename, particles, N);

    for (int ts = 0; ts < nsteps; ts++) {
        for (int i = 0; i < N; i++) {
            Vector total_force = {0.0, 0.0};
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    Vector f_ij = force_particle(&particles[i], &particles[j], N);
                    total_force.x += f_ij.x;
                    total_force.y += f_ij.y;
                }
            }
            if (particles[i].mass != 0) {
                particles[i].acceleration.x = total_force.x / particles[i].mass;
                particles[i].acceleration.y = total_force.y / particles[i].mass;
            } else {
                particles[i].acceleration.x = 0;
                particles[i].acceleration.y = 0;
            }
            particles[i].velocity.x += delta_t * particles[i].acceleration.x;
            particles[i].velocity.y += delta_t * particles[i].acceleration.y;
            particles[i].position.x += delta_t * particles[i].velocity.x;
            particles[i].position.y += delta_t * particles[i].velocity.y;
        }
    }

    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        perror("Error opening the file!\n");
        free(particles);
        fclose(output_file);
        return 1;
    }

    // Calculate the total size of data to write in bytes
    size_t data_size = N * sizeof(double) * 6;

    // Allocate memory for the data
    double* data_buffer = malloc(data_size);
    if (data_buffer == NULL) {
        perror("Memory allocation failed.\n");
        free(particles);
        fclose(output_file);
        return 1;
    }

    // Copy the data to the buffer
    for (int i = 0; i < N; i++) {
        memcpy(data_buffer + i * 6, &particles[i].position.x, sizeof(double));
        memcpy(data_buffer + i * 6 + 1, &particles[i].position.y, sizeof(double));
        memcpy(data_buffer + i * 6 + 2, &particles[i].mass, sizeof(double));
        memcpy(data_buffer + i * 6 + 3, &particles[i].velocity.x, sizeof(double));
        memcpy(data_buffer + i * 6 + 4, &particles[i].velocity.y, sizeof(double));
        memcpy(data_buffer + i * 6 + 5, &particles[i].brightness, sizeof(double));
    }

    // Write the data buffer to the file
    fwrite(data_buffer, sizeof(double), N * 6, output_file);

    // Free the allocated memory
    free(data_buffer);
    free(particles);
    fclose(output_file);

    return 0;
}

