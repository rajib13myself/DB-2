#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NX 100  // Number of grid points in x-direction
#define NY 100  // Number of grid points in y-direction
#define NSTEPS 1000  // Number of time steps
#define DT 0.001  // Time step
#define DX 0.1  // Grid spacing in x-direction
#define DY 0.1  // Grid spacing in y-direction
#define ALPHA 0.1  // Diffusion coefficient

// Function to initialize the temperature distribution
void initialize(double** u, int nx, int ny) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            u[i][j] = 0.0;  // Initial temperature
        }
    }
    
    // Boundary conditions
    for (int i = 0; i < nx; i++) {
        u[i][0] = 100.0;  // Bottom boundary
        u[i][ny - 1] = 0.0;  // Top boundary
    }
    for (int j = 0; j < ny; j++) {
        u[0][j] = 75.0;  // Left boundary
        u[nx - 1][j] = 50.0;  // Right boundary
    }
}

// Function to perform one time step using finite difference method
void step(double** u, double** un, int nx, int ny, double dt, double dx, double dy, double alpha) {
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            un[i][j] = u[i][j] + alpha * dt * (
                (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / (dx * dx) +
                (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / (dy * dy)
            );
           // printf("u:%lf", u)
        }
        printf("u:%lf", u);
    }

    // Update u for the next time step
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            u[i][j] = un[i][j];
        }
    }
    printf("After update u: %lf", u );
}

int main() {
    // Allocate memory for temperature grid
    double **u = (double **)malloc(NX * sizeof(double *));
    double **un = (double **)malloc(NX * sizeof(double *));
    for (int i = 0; i < NX; i++) {
        u[i] = (double *)malloc(NY * sizeof(double));
        un[i] = (double *)malloc(NY * sizeof(double));
    }

    // Initialize temperature distribution
    initialize(u, NX, NY);

    // Perform time integration
    for (int t = 0; t < NSTEPS; t++) {
        step(u, un, NX, NY, DT, DX, DY, ALPHA);
    }

    // Output results or perform further analysis

    // Free memory
    for (int i = 0; i < NX; i++) {
        free(u[i]);
        free(un[i]);
    }
    free(u);
    free(un);

    return 0;
}

