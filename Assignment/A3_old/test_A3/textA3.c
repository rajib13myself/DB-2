#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// Define particle structure
typedef struct {
    double x, y;    // Position
    double velx, vely;  // Velocity
    double mass;    // Mass
    double brightness;  // Brightness
} Particle;

// Function to read initial configuration from file
Particle* read_initial_configuration(const char* filename, int N);

// Function to simulate particle motion for given number of timesteps
void simulate_motion(Particle* particles, int N, int nsteps, double delta_t, bool graphics);

// Function to save final positions and velocities to result file
void save_result(const char* filename, Particle* particles, int N);

int main(int argc, char *argv[]) {
    // Check if correct number of arguments is provided
    if (argc != 6) {
        printf("Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
        return 1;
    }

    // Parse input arguments
    int N = atoi(argv[1]);
    const char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    bool graphics = atoi(argv[5]);

    // Read initial configuration from file
    Particle* particles = read_initial_configuration(filename, N);
    if (particles == NULL) {
        printf("Error reading initial configuration from file.\n");
        return 1;
    }

    // Simulate particle motion
    simulate_motion(particles, N, nsteps, delta_t, graphics);

    // Save final positions and velocities to result file
    save_result("result.gal", particles, N);

    // Free allocated memory
    free(particles);

    return 0;
}

Particle* read_initial_configuration(const char* filename, int N) {
    // Open file for reading
    
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        printf("Error opening file: %s\n", filename);
        return NULL;
    }

    // Allocate memory for particles
    Particle* particles = malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Memory allocation failed.\n");
        fclose(file);
        return NULL;
    }

    // Read particle data from file
    size_t num_read = fread(particles, sizeof(Particle), N, file);
    fclose(file);

    // Check if correct number of particles read
    if (num_read != N) {
        printf("Error reading particle data from file.\n");
        free(particles);
        return NULL;
    }

    return particles;
}

void simulate_motion(Particle* particles, int N, int nsteps, double delta_t, bool graphics) {
    // Simulation logic goes here
    // This is where you would implement the simulation of particle motion
}

void save_result(const char* filename, Particle* particles, int N) {
    // Open file for writing
    FILE* file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file: %s\n", filename);
        return;
    }

    // Write particle data to file
    fwrite(particles, sizeof(Particle), N, file);

    fclose(file);
}

