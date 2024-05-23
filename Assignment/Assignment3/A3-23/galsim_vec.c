#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void matvec_ref(double **mat_a, double *vec_b, double *vec_c, int SIZE) {
    int i, j;

    for (i = 0; i < SIZE; i++) {
        double d = 0.0;
        for (j = 0; j < SIZE; j++) {
            d += mat_a[i][j] * vec_b[j];
        }
        vec_c[i] = d;
    }
}

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

    // Allocate memory for particle properties
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *mass = malloc(N * sizeof(double));
    double *vx = malloc(N * sizeof(double));
    double *vy = malloc(N * sizeof(double));
    double *brightness = malloc(N * sizeof(double));

    if (x == NULL || y == NULL || mass == NULL || vx == NULL || vy == NULL || brightness == NULL) {
        printf("Error: Memory allocation failed.\n");
        return 0;
    }

    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
        free(x);
        free(y);
        free(mass);
        free(vx);
        free(vy);
        free(brightness);
        return 0;
    }

    // Read particle data
    for (int i = 0; i < N; i++) {
        fread(&x[i], sizeof(double), 1, input_file);
        fread(&y[i], sizeof(double), 1, input_file);
        fread(&mass[i], sizeof(double), 1, input_file);
        fread(&vx[i], sizeof(double), 1, input_file);
        fread(&vy[i], sizeof(double), 1, input_file);
        fread(&brightness[i], sizeof(double), 1, input_file);
    }

    fclose(input_file);

    // Simulation loop
    for (int step = 0; step < nsteps; step++) {
        // Calculate forces
        double *forces_x = calloc(N, sizeof(double));
        double *forces_y = calloc(N, sizeof(double));
        double *acc_x = calloc(N, sizeof(double));
        double *acc_y = calloc(N, sizeof(double));

        if (forces_x == NULL || forces_y == NULL || acc_x == NULL || acc_y == NULL) {
            printf("Error: Memory allocation failed.\n");
            free(forces_x);
            free(forces_y);
            free(acc_x);
            free(acc_y);
            return 0;
        }

        double **mat_a = malloc(N * sizeof(double *));
        for (int i = 0; i < N; i++) {
            mat_a[i] = calloc(N, sizeof(double));
            if (mat_a[i] == NULL) {
                printf("Error: Memory allocation failed.\n");
                free(forces_x);
                free(forces_y);
                free(acc_x);
                free(acc_y);
                return 0;
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    double rx = x[i] - x[j];
                    double ry = y[i] - y[j];
                    double r = sqrt(rx * rx + ry * ry) + 1e-3; // Add epsilon to avoid division by zero
                    double f_ij = -(100.0 / N) * mass[i] * mass[j] / (r * r * r);
                    mat_a[i][j] = f_ij / mass[i];
                }
            }
        }

        matvec_ref(mat_a, vx, acc_x, N);
        matvec_ref(mat_a, vy, acc_y, N);

        // Update velocities
        for (int i = 0; i < N; i++) {
            vx[i] += acc_x[i] * delta_t;
            vy[i] += acc_y[i] * delta_t;
        }

        free(forces_x);
        free(forces_y);
        free(acc_x);
        free(acc_y);

        // Update positions
        for (int i = 0; i < N; i++) {
            x[i] += vx[i] * delta_t;
            y[i] += vy[i] * delta_t;
        }

        for (int i = 0; i < N; i++) {
            free(mat_a[i]);
        }
        free(mat_a);
    }

    // Write results to file
    FILE *output_file = fopen("result.gal", "wb");
    if (output_file == NULL) {
        printf("Error: Unable to open output file.\n");
        free(x);
        free(y);
        free(mass);
        free(vx);
        free(vy);
        free(brightness);
        return 0;
    }

    // Write particle data
    for (int i = 0; i < N; i++) {
        fwrite(&x[i], sizeof(double), 1, output_file);
        fwrite(&y[i], sizeof(double), 1, output_file);
        fwrite(&mass[i], sizeof(double), 1, output_file);
        fwrite(&vx[i], sizeof(double), 1, output_file);
        fwrite(&vy[i], sizeof(double), 1, output_file);
        fwrite(&brightness[i], sizeof(double), 1, output_file);
    }

    fclose(output_file);

    // Free allocated memory
    free(x);
    free(y);
    free(mass);
    free(vx);
    free(vy);
    free(brightness);

    return 0;
}
