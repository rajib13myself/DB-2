#include <stdio.h>
#include <omp.h>

#define PLATE_LENGTH 50
#define MAX_ITER_TIME 1000

#define ALPHA  2.0
#define DELTA_X 1.0

// Calculated parameters
#define DELTA_T ((DELTA_X * DELTA_X) / (4 * ALPHA))
#define GAMMA ((ALPHA * DELTA_T) / (DELTA_X * DELTA_X))

// Initialize solution: the grid of u(k, i, j)
double u[MAX_ITER_TIME][PLATE_LENGTH][PLATE_LENGTH];

// Function to set initial conditions
void set_InitialConditions() {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < MAX_ITER_TIME; i++) {
        for (int j = 0; j < PLATE_LENGTH; j++) {
            for (int k = 0; k < PLATE_LENGTH; k++) {
                u[i][j][k] = 0.0;
            }
        }
    }
}

// Function to set boundary conditions
void set_BoundaryConditions() {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < MAX_ITER_TIME; i++) {
        for (int j = 0; j < PLATE_LENGTH; j++) {
            u[i][j][0] = 0.0;  // left boundary
            u[i][0][j] = 0.0;  // bottom boundary
            u[i][PLATE_LENGTH - 1][j] = 100.0;  // top boundary
            u[i][j][PLATE_LENGTH - 1] = 0.0;  // right boundary
        }
    }
}

// Function to perform the explicit finite difference method
void solve_HeatEquation() {
    #pragma omp parallel for collapse(3)
    for (int k = 1; k < MAX_ITER_TIME; k++) {
        for (int i = 1; i < PLATE_LENGTH - 1; i++) {
            for (int j = 1; j < PLATE_LENGTH - 1; j++) {
                u[k][i][j] = u[k - 1][i][j] + GAMMA * (u[k - 1][i + 1][j] + u[k - 1][i - 1][j] + u[k - 1][i][j + 1] + u[k - 1][i][j - 1] - 4 * u[k - 1][i][j]);
            }
        }
    }
}

int main() {
    set_InitialConditions();
    set_BoundaryConditions();
    solve_HeatEquation();

    // Save the final state of the plate to a text file
    FILE *file = fopen("result.txt", "w");
    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }
    fprintf(file, "Final state of the plate:\n");
    for (int i = 0; i < PLATE_LENGTH; i++) {
        for (int j = 0; j < PLATE_LENGTH; j++) {
            fprintf(file, "%.2f ", u[MAX_ITER_TIME - 1][i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);

    return 0;
}
