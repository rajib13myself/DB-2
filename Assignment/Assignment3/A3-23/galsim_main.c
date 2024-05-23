#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALIGNMENT 32 // Adjust alignment for aligned memory allocation

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    double brightness; // brightness
} Particle;

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
        return 0;
    }

    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);

    // Allocate memory for particles
    Particle *particles = (Particle*)aligned_alloc(ALIGNMENT, N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 0;
    }

    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        fclose(input_file);
        free(particles);
        return 0;
    }

    for (int i = 0; i < N; i++) {
        fread(&particles[i].x, sizeof(double), 1, input_file);
        fread(&particles[i].y, sizeof(double), 1, input_file);
        fread(&particles[i].mass, sizeof(double), 1, input_file);
        fread(&particles[i].vx, sizeof(double), 1, input_file);
        fread(&particles[i].vy, sizeof(double), 1, input_file);
        fread(&particles[i].brightness, sizeof(double), 1, input_file);
    }

    fclose(input_file);

    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        register double ax, ay;
        for (int i = 0; i < N; i++) {
            ax = 0.0;  ay = 0.0;
            for (int j = 0; j < N; j+=2) { // Unroll by 2
                register double rx1, ry1, rx2, ry2, r1, r2;
                register double f_ij1, f_ij2;
                
                rx1 = (particles[j].x - particles[i].x);
                ry1 = (particles[j].y - particles[i].y);
                r1 = sqrt(rx1 * rx1 + ry1 * ry1) + 1e-3;
                if (j != i) {
                    f_ij1 = (100.0/N) * (particles[i].mass) * particles[j].mass /(r1 * r1 * r1);
                    ax += f_ij1 * rx1 / particles[i].mass;
                    ay += f_ij1 * ry1 / particles[i].mass;
                }

                rx2 = (particles[j+1].x - particles[i].x);
                ry2 = (particles[j+1].y - particles[i].y);
                r2 = sqrt(rx2 * rx2 + ry2 * ry2) + 1e-3;
                if (j+1 != i) {
                    f_ij2 = (100.0/N) * (particles[i].mass) * particles[j+1].mass /(r2 * r2 * r2);
                    ax += f_ij2 * rx2 / particles[i].mass;
                    ay += f_ij2 * ry2 / particles[i].mass;
                }
            }
            particles[i].vx += ax * delta_t;
            particles[i].vy += ay * delta_t;
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
        fclose(output_file);
        free(particles);
        return 0;
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

    return 0;
}
