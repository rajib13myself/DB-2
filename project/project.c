#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>

#define TGX 50  // Number of grid points in x-direction
#define TGY 50  // Number of grid points in y-direction
#define TSTEPS 10000 // Number of time steps
#define TS 0.001  // Time step
#define GSX 0.1  // Grid spacing in x-direction
#define GSY 0.1  // Grid spacing in y-direction
#define ALPHA 0.1  // Diffusion coefficient

// Precomputed constants
const double GSX_GSY = GSX * GSY;
const double ALPHA_TS = ALPHA * TS;
const double GSX2 = GSX * GSX;
const double GSY2 = GSY * GSY;

// Function to allocate memory for the temperature grid
double* allocate_temp_grid(int size) {
    double* grid = (double*)calloc(size, sizeof(double));
    if (grid == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    return grid;
}

// Function to free memory allocated for the temperature grid
void free_temp_grid(double* grid) {
    free(grid);
}

// Function to initialize the temperature distribution
void initialize_pde(double* u) {
    // Initialize grid with 0s
    for (int i = 0; i < TGX * TGY; i++) {
        u[i] = 0.0;  // Initial temperature
    }

    // Initialize Boundary conditions
    for (int i = 0; i < TGX; i++) {
        u[i * TGY] = 0.0;          // Bottom boundary
        u[i * TGY + TGY - 1] = 50.0; // Top boundary
    }
    for (int j = 0; j < TGY; j++) {
        u[j] = 0.0;               // Left boundary
        u[(TGX - 1) * TGY + j] = 50.0; // Right boundary
    }
}

// Function to perform one time step using finite difference method
void timestep_fd(double** u, double** un, FILE* file) {
    #pragma omp parallel for
    for (int i = 1; i < TGX - 1; i++) {
        #pragma omp simd
        for (int j = 1; j < TGY - 1; j++) {
            (*un)[i * TGY + j] = (*u)[i * TGY + j] + ALPHA_TS * (
                ((*u)[(i+1) * TGY + j] - 2.0 * (*u)[i * TGY + j] + (*u)[(i-1) * TGY + j]) / GSX2 +
                ((*u)[i * TGY + (j+1)] - 2.0 * (*u)[i * TGY + j] + (*u)[i * TGY + (j-1)]) / GSY2
            );
        }
    }

    // Swap u and un for the next time step
    double* temp = *u;
    *u = *un;
    *un = temp;

    // Write temperature data to file
    for (int i = 0; i < TGX; i++) {
        for (int j = 0; j < TGY; j++) {
            fprintf(file, "%lf ", (*u)[i * TGY + j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
    fflush(file); // Flush the file stream to ensure data is written to the file immediately
}

int main() {
    // Allocate memory for the temperature grids
    double* u = allocate_temp_grid(TGX * TGY);
    double* un = allocate_temp_grid(TGX * TGY);

    // Initialize temperature distribution
    initialize_pde(u);

    // Open file for writing temperature data
    FILE *file = fopen("temperature_data.txt", "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    // Perform time integration
    for (int t = 0; t < TSTEPS; t++) {
        timestep_fd(&u, &un, file);
    }

    fclose(file);

    // Free allocated memory
    free_temp_grid(u);
    free_temp_grid(un);

    return 0;
}
