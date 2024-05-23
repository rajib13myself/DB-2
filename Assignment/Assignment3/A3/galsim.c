#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
        fclose(input_file);
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
        double rx, ry, r, f_ij, f_ji;
        double ax, ay;
        //double *forces_x = calloc(N, sizeof(double));
        //double *forces_y = calloc(N, sizeof(double));
        /*
        if (forces_x == NULL || forces_y == NULL) {
            printf("Error: Memory allocation failed.\n");
            free(forces_x);
            free(forces_y);
            return 0;
        }*/

        // Calculate forces
        for (int i = 0; i < N; i++) {
            //forces_x[i] = 0.0; forces_y[i] = 0.0;
            ax = 0.0; ay = 0.0;
            f_ij = (100.0 / N) * mass[i];
            for (int j = 0; j < N; j++) {  
                rx = x[i] - x[j];
                ry = y[i] - y[j];
                r = sqrt(rx * rx + ry * ry) + 1e-3; // Add epsilon to avoid division by zero
                if(j!=i){
                    //f_ij = (100.0 / N) * mass[i] * mass[j] / (r * r * r);
                    f_ji = -f_ij * mass[j] / (r * r * r);
                    //f_ji = -f_ij; // Newton's third law
                    ax += f_ji * rx / mass[i];
                    ay += f_ji * ry / mass[i];
                    //printf("\nax and ay: %lf %lf", ax, ay);
                }
                
            }

            // Update velocities
            vx[i] += ax * delta_t;
            vy[i] += ay * delta_t;
        }
        
        // Free memory for forces
        //free(forces_x);
        //free(forces_y);

        // Update positions
        for (int i = 0; i < N; i++) {
            x[i] += vx[i] * delta_t;
            y[i] += vy[i] * delta_t;
        }
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
        fclose(output_file);
        return 0;
    }

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
