#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <emmintrin.h> // Include SIMD intrinsics header for SSE

#define ALIGNMENT 16 // Adjust alignment for aligned memory allocation

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
        return 1;
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
        double ax, ay;
        for (int i = 0; i < N; i++) {
            ax = 0.0;  ay = 0.0;
            for (int j = 0; j < N; j+=2) { // Vectorize by 2 (example)
                __m128d rx, ry, r_squared, f_ij, ax_part, ay_part;

                // Load particle attributes from memory into SIMD registers
                rx = _mm_load_pd(&particles[j].x); // Load 2 x-coordinates
                ry = _mm_load_pd(&particles[j].y); // Load 2 y-coordinates

                // Compute distance squared (rx * rx + ry * ry)
                r_squared = _mm_add_pd(_mm_mul_pd(rx, rx), _mm_mul_pd(ry, ry));

                // Compute inverse of distance cubed (1 / r^3)
                f_ij = _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(r_squared));
                f_ij = _mm_mul_pd(f_ij, _mm_mul_pd(f_ij, f_ij));

                // Compute forces (f_ij = (100.0/N) * particles[i].mass * particles[j].mass * f_ij)
                f_ij = _mm_mul_pd(_mm_set1_pd(100.0 / N), f_ij);
                f_ij = _mm_mul_pd(f_ij, _mm_mul_pd(_mm_set1_pd(particles[i].mass), _mm_set_pd(particles[j].mass[j+1], particles[j].mass[j])));

                // Compute force components along x and y directions
                ax_part = _mm_mul_pd(f_ij, rx);
                ay_part = _mm_mul_pd(f_ij, ry);

                // Accumulate forces along x and y directions
                ax_part = _mm_hadd_pd(ax_part, ax_part);
                ay_part = _mm_hadd_pd(ay_part, ay_part);

                // Store accumulated forces to scalar variables
                double ax_vals[2], ay_vals[2];
                _mm_store_pd(ax_vals, ax_part);
                _mm_store_pd(ay_vals, ay_part);
                ax += ax_vals[0] + ax_vals[1];
                ay += ay_vals[0] + ay_vals[1];
            }

            // Update velocities
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
