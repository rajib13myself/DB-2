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

// Function to allocate memory for the temperature grid
double** allocate_temp_grid(int rows, int cols) {
    double** grid = (double**)malloc(rows * sizeof(double*));
    if (grid == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; i++) {
        grid[i] = (double*)malloc(cols * sizeof(double));
        if (grid[i] == NULL) {
            perror("Memory allocation failed");
            exit(EXIT_FAILURE);
        }
    }
    return grid;
}

// Function to free memory allocated for the temperature grid
void free_temp_grid(double** grid, int rows) {
    for (int i = 0; i < rows; i++) {
        free(grid[i]);
    }
    free(grid);
}

// Function to initialize the temperature distribution
void initialize_pde(double** u) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < TGX; i++) {
        for (int j = 0; j < TGY; j++) {
            u[i][j] = 0.0;  // Initial temperature
        }
    }

    // Initialize Boundary conditions
    #pragma omp parallel for
    for (int i = 0; i < TGX; i++) {
        u[i][0] = 50.0;    // Bottom boundary
        u[i][TGY - 1] = 0.0; // Top boundary
    }
    #pragma omp parallel for
    for (int j = 0; j < TGY; j++) {
        u[0][j] = 50.0;     // Left boundary
        u[TGX - 1][j] = 50.0; // Right boundary
    }
}

// Function to perform one time step using finite difference method
void timestep_fd(double** u, double** un) {
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < TGX - 1; i++) {
        for (int j = 1; j < TGY - 1; j++) {
            un[i][j] = u[i][j] + ALPHA * TS * (
                (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / (GSX * GSX) +
                (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / (GSY * GSY)
            );
        }
    }

    // Update u for the next time step
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < TGX - 1; i++) {
        for (int j = 1; j < TGY - 1; j++) {
            u[i][j] = un[i][j];
        }
    }
}

int main() {
    // Allocate memory for the temperature grids
    double** u = allocate_temp_grid(TGX, TGY);
    double** un = allocate_temp_grid(TGX, TGY);

    // Initialize temperature distribution
    initialize_pde(u);
    
   //Open a file for write the output result of time integration
   char filename[50]; 	
   for( int t = 0; t < TSTEPS; t++) {
   	timestep_fd(u, un);
	sprintf(filename, "temperature_%d.txt", t);
	FILE *output_file = fopen(filename, "w");
	if(output_file == NULL) {
		perror("Error in opening file");
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < TGX; i++) {
		for(int j = 0; j < TGY; j++) {
			fprintf(output_file, "%lf ", u[i][j]);
		}
		fprintf(output_file, "\n");
	}
	fclose(output_file);
    }


    // Free allocated memory
    free_temp_grid(u, TGX);
    free_temp_grid(un, TGX);

    //pclose(gplot);
    return 0;
}

