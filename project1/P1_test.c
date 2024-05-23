#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>

#define NGX 50  // Number of grid points in x-direction
#define NGY 50  // Number of grid points in y-direction
#define NTSTEPS 10000 // Number of time steps
#define DTS 0.001  // Time step
#define DGX 0.1  // Grid spacing in x-direction
#define DGY 0.1  // Grid spacing in y-direction
#define ALPHA 0.1  // Diffusion coefficient

// Function to initialize the temperature distribution
void initialize_pde(double u[NGX][NGY]) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NGX; i++) {
        for (int j = 0; j < NGY; j++) {
            u[i][j] = 0.0;  // Initial temperature
        }
    }

    // Initialize Boundary conditions
    #pragma omp parallel for
    for (int i = 0; i < NGX; i++) {
        u[i][0] = 50.0;    // Bottom boundary
        u[i][NGY - 1] = 0.0; // Top boundary
    }
    #pragma omp parallel for
    for (int j = 0; j < NGY; j++) {
        u[0][j] = 50.0;     // Left boundary
        u[NGX - 1][j] = 50.0; // Right boundary
    }
}

// Function to perform one time step using finite difference method
void step_fd(double u[NGX][NGY], double un[NGX][NGY]) {
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NGX - 1; i++) {
        for (int j = 1; j < NGY - 1; j++) {
            un[i][j] = u[i][j] + ALPHA * DTS * (
                (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / (DGX * DGX) +
                (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / (DGY * DGY)
            );
        }
    }

    // Update u for the next time step
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NGX - 1; i++) {
        for (int j = 1; j < NGY - 1; j++) {
            u[i][j] = un[i][j];
        }
    }
}

int main() {
    double u[NGX][NGY]; // Temperature grid
    double un[NGX][NGY]; // Temporary grid for next time step

    // Initialize temperature distribution
    initialize_pde(u);



    // Perform time integration
    FILE *gplot = popen("gnuplot -persist", "w");
    fprintf(gplot, "set pm3d map\n");
    fprintf(gplot, "set xlabel 'X'\n");
    fprintf(gplot, "set ylabel 'Y'\n");
    fprintf(gplot, "set zlabel 'Temperature'\n");
    fprintf(gplot, "set xrange [0:%d]\n", NGX - 1);
    fprintf(gplot, "set yrange [0:%d]\n", NGY - 1);
    
  //  #pragma omp parallel
   //{   
    for (int t = 0; t < NTSTEPS; t++) {
        step_fd(u, un);
        fprintf(gplot, "splot '-' matrix with image\n");
        for (int i = 0; i < NGX; i++) {
            for (int j = 0; j < NGY; j++) {
                fprintf(gplot, "%lf ", u[i][j]);
            }
            fprintf(gplot, "\n");
        }
        fprintf(gplot, "e\n");
        fflush(gplot);
        usleep(10000); // Sleep for visualization
    }
   //}
    
    pclose(gplot);
    return 0;
}
