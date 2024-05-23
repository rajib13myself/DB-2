#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double data[6]; // 1D array for position (x, y), velocity (vx, vy), mass, brightness
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
    Particle *particles = (Particle*)malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 0;
    }

    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        free(particles);
        return 0;
    }

    for (int i = 0; i < N ; i++) {
        fread(particles[i].data, sizeof(double), 6, input_file);
        // Print values for debugging
        //printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].data[0], particles[i].data[1], particles[i].data[2], particles[i].data[3], particles[i].data[4], particles[i].data[5]);
    }

    fclose(input_file);

    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        /* //Allocate memory for forces
        double *forces_x = (double*)malloc(N * sizeof(double));
        double *forces_y = (double*)malloc(N * sizeof(double));
        if (forces_x == NULL || forces_y == NULL) {
            printf("Error: Memory allocation failed.\n");
            free(forces_x);
            free(forces_y);
            return 0;
        }*/
        register double r, rx, ry, f_ij;
        double forces_x, forces_y;
        // Zero-initialize forces arrays
        for (int i = 0; i < N; i++) {
            forces_x = 0.0;
            forces_y = 0.0;
            //double rx= 0.0, ry = 0.0;
            //double ax_ij[i] = 0.0;
            for (int j = 0; j < N; j++) {
                    rx = particles[i].data[0] - particles[j].data[0];
                    ry = particles[i].data[1] - particles[j].data[1];
                    //Debug value of x, y
                    //printf("rx and ry for i,j: %lf , %lf\n", rx, ry);
                    r = sqrt(rx * rx + ry * ry) + 1e-3; // Add epsilon to avoid division by zero
                    if (j != i) {
                        f_ij = -(100.0 / N) * particles[i].data[2] * particles[j].data[2] / (r * r * r);
                        forces_x += f_ij * rx / particles[i].data[2];
                        forces_y += f_ij * ry / particles[i].data[2];
                }
            }
            // Calculate acceleration and update velocities
            
            particles[i].data[3] += forces_x * delta_t;
            particles[i].data[4] += forces_y * delta_t;
            //printf("%f  %f\n", forces_x, forces_y);
            
        }

        // Free memory for forces
        //free(forces_x);
        //free(forces_y);

        // Update positions
        for (int i = 0; i < N; i++) {
            particles[i].data[0] += particles[i].data[3] * delta_t;
            particles[i].data[1] += particles[i].data[4] * delta_t;
        } 


    }

    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        free(particles);
        return 0;
    }

    for (int i = 0; i < N; i++) {
        fwrite(particles[i].data, sizeof(double), 6, output_file);
        // Print values for debugging
        //printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].data[0], particles[i].data[1], particles[i].data[2], particles[i].data[3], particles[i].data[4], particles[i].data[5]);
    }

    fclose(output_file);

    // Free allocated memory
    free(particles);

    return 0;
}
