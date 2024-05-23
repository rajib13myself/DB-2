#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define G 100.0    // Updated gravitational constant as per the instructions
#define EPSILON 1e-3
#define DT 1e-5

typedef struct {
    double x, y;
} Vector;

typedef struct {
    double mass;
    Vector position;
    Vector velocity;
    Vector acceleration;
    double brightness;
} Particle;

Vector vector_cal(Particle p1, Particle p2) {
    Vector r1_p1p2;
    r1_p1p2.x = p2.position.x - p1.position.x;
    r1_p1p2.y = p2.position.y - p1.position.y;
    return r1_p1p2;
}

double distance_cal(Particle p1, Particle p2) {
    double dx = p2.position.x - p1.position.x;
    double dy = p2.position.y - p1.position.y;
    double dist_val = sqrt(dx * dx + dy * dy);
    if (dist_val <= 0) {
        return 0.0001;
    }
    return dist_val;
}

Vector distance_normalize(Vector r_ij) {
    double dist_fn = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y);
    Vector r2_p1p2 = {0};
    if (dist_fn != 0) {
        r2_p1p2.x = r_ij.x / dist_fn;
        r2_p1p2.y = r_ij.y / dist_fn;
    }
    return r2_p1p2;
}

Particle* allocate_particles(int N) {
    Particle* particles = malloc(N * sizeof(Particle));
    if (particles == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    return particles;
}

void read_particles(const char* filename, Particle* particles, int N) {
    FILE* read_file = fopen(filename, "rb");
    if (read_file == NULL) {
        fprintf(stderr, "Error opening the file.\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fread(&particles[i], sizeof(Particle), 1, read_file);
    }
    fclose(read_file);
}

Vector force_particle(Particle *p1, Particle *p2) {
    Vector r_ij = vector_cal(*p1, *p2);
    double dist_fn = distance_cal(*p1, *p2);
    if (dist_fn == 0 || p1->mass == 0 || p2->mass == 0) {
        Vector def_force = {0, 0};
        return def_force;
    }
    int N;
    Vector r_hat = distance_normalize(r_ij);
    double f = -(100.0 / N) * ((p1->mass * p2->mass) / ((dist_fn + EPSILON) * (dist_fn + EPSILON) * (dist_fn + EPSILON)));
    Vector f_p1p2 = {f * r_hat.x, f * r_hat.y};
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
        // Calculate forces on each particle
        for (int i = 0; i < N; i++) {
            Vector total_force = {0.0, 0.0};
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    Vector f_ij = force_particle(&particles[i], &particles[j]);
                    total_force.x += f_ij.x;
                    total_force.y  += f_ij.y; // Add force components
                }
            }
            // Update Acceleration
            if (particles[i].mass != 0) {
                particles[i].acceleration.x = total_force.x / particles[i].mass;
                particles[i].acceleration.y = total_force.y / particles[i].mass;
            } else {
                particles[i].acceleration.x = 0;
                particles[i].acceleration.y = 0;
            }

            // Update Velocity using symplectic Euler method
            particles[i].velocity.x += delta_t * particles[i].acceleration.x;
            particles[i].velocity.y += delta_t * particles[i].acceleration.y;

            // Update Position using symplectic Euler method
            particles[i].position.x += delta_t * particles[i].velocity.x;
            particles[i].position.y += delta_t * particles[i].velocity.y;
        }
    }

    // Output results to a file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        perror("Error opening the file!\n");
        free(particles);
        fclose(output_file);
        return 1;
    }

    /*
    for (int i = 0; i < N; i++) {
        fwrite(&particles[i], sizeof(Particle), 1, output_file);
    }

    fclose(output_file);
    free(particles);
    return 0;
*/

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

