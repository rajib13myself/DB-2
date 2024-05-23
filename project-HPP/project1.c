#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#define TGX 50  // Number of grid points in x-direction
#define TGY 50  // Number of grid points in y-direction
#define TSTEPS 10000 // Number of time steps
#define TS 0.001  // Time step
#define GSX 0.1  // Grid spacing in x-direction
#define GSY 0.1  // Grid spacing in y-direction
#define ALPHA 0.1  // Diffusion coefficient
#define UPDATE_INTERVAL 100  // Update plot every UPDATE_INTERVAL time steps

// Precomputed constants
const double GSX_GSY = GSX * GSY;
const double ALPHA_TS_GSX2 = ALPHA * TS / (GSX * GSX);
const double ALPHA_TS_GSY2 = ALPHA * TS / (GSY * GSY);

// Function to allocate memory for the temperature grid
double* allocate_temp_grid(int rows, int cols) {
    double* grid = (double*)calloc(rows * cols, sizeof(double));
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

// Function to perform one time step using finite difference method
void timestep_fd(double* u, double* un) {
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < TGX - 1; i++) {
        #pragma omp simd
        for (int j = 1; j < TGY - 1; j++) {
            un[i * TGY + j] = u[i * TGY + j] + ALPHA_TS_GSX2 * (
                (u[(i + 1) * TGY + j] - 2.0 * u[i * TGY + j] + u[(i - 1) * TGY + j]) +
                ALPHA_TS_GSY2 * (u[i * TGY + (j + 1)] - 2.0 * u[i * TGY + j] + u[i * TGY + (j - 1)])
            );
        }
    }
    // Swap u and un for the next time step
    double* temp = u;
    u = un;
    un = temp;
}

int main(int argc, char *argv[]) {
    int graphics_enabled = 1; // By default, graphics is enabled
    if (argc == 2 && atoi(argv[1]) == 0) {
        graphics_enabled = 0;
    }

    // Allocate memory for the temperature grids
    double* u = allocate_temp_grid(TGX, TGY);
    double* un = allocate_temp_grid(TGX, TGY);

    // Initialize temperature distribution
    // Initialization is not needed as calloc initializes memory with zeros

    // Perform time integration
    FILE *file = fopen("temperature_data.txt", "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    for (int t = 0; t < TSTEPS; t++) {
        timestep_fd(u, un);
        // Write temperature data to file after every UPDATE_INTERVAL time steps
        if (t % UPDATE_INTERVAL == 0) {
            for (int i = 0; i < TGX; i++) {
                for (int j = 0; j < TGY; j++) {
                    fprintf(file, "%lf ", u[i * TGY + j]);
                }
                fprintf(file, "\n");
            }
            fprintf(file, "\n");
            fflush(file); // Flush the file stream to ensure data is written to the file immediately
        }
    }

    fclose(file);

    // Plot data from file using gnuplot
    if (graphics_enabled) {
        FILE *gplot = popen("gnuplot -persist", "w");
        fprintf(gplot, "set pm3d map\n");
        fprintf(gplot, "set xlabel 'X'\n");
        fprintf(gplot, "set ylabel 'Y'\n");
        fprintf(gplot, "set zlabel 'Temperature'\n");
        fprintf(gplot, "set xrange [0:%d]\n", TGX - 1);
        fprintf(gplot, "set yrange [0:%d]\n", TGY - 1);
        fprintf(gplot, "splot 'temperature_data.txt' matrix with image\n");
        pclose(gplot);
    }

    // Free allocated memory
    free_temp_grid(u);
    free_temp_grid(un);

    return 0;
}
