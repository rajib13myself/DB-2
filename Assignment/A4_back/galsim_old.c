#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>

#define G 100.0 // gravitational constant adjusted by number of particles
//#define EPSILON 1e-3 // small number for stability

//Define mutex for thread syncronization
pthread_mutex_t mutex;

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double mass; // mass
    //    double ax, ay; //acceleration
    double brightness; // brightness
} Particle;

typedef struct {
    int start;
    int end;
    int N;
    Particle *particles;
} ThreadArgs;

/*
// Function to calculate distance between two particles
double distance(Particle *p1, Particle *p2) {
    return sqrt(((p2->x - p1->x) * (p2->x - p1->x))  + ((p2->y - p1->y) * (p2->y - p1->y)));

}


// Function to calculate total force acting on a particle
double force(Particle *p1, Particle *p2,int N) {
    double r = distance(p1, p2);
    double r_e = r + EPSILON;
    return -(100.0/N) * p1->mass * p2->mass /(r_e * r_e * r_e);
}*/

//Function to be executed by each Thread
void* simulation_force(void* args) {
    ThreadArgs* threadArgs = (ThreadArgs*)args;
    Particle* particles = threadArgs->particles;
    int start = threadArgs->start;
    int end = threadArgs->end;
    int N = threadArgs->N;
    double delta_t = *((double*)args);
    register double f_ij;
    register double ax, ay;
    register double r, rx, ry;
    for (int i = start; i < end; i ++) {
       ax = 0.0; ay = 0.0;
       for (int j = 0; j < N; j++) {
          rx = particles[j].x - particles[i].x;
          ry = particles[j].y - particles[i].y;
	  r = sqrt(rx * rx + ry * ry) + 1e-3;
          if ( j != i) {
             //f_ij = force(&particles[i], &particles[j],N);
	     f_ij = (100.0/N) * particles[i].mass * particles[j].mass /(r * r * r);
             ax += f_ij * rx/particles[i].mass;
             ay += f_ij * ry/particles[i].mass;
          }
                
       }
       particles[i].vx += ax * delta_t;
       particles[i].vy += ay * delta_t;
       // printf("%f   %f\n",ax,ay);
       // double dummy;
       //scanf("stop:",&dummy);
    }
    return NULL;
} 

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Usage: %s N filename nsteps delta_t graphics n_thread\n", argv[0]);
        return 1;
    }
    
    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    int num_threads = atoi(argv[6]);
    
    //Initialize mutex
    pthread_mutex_init(&mutex, NULL);

    
    // Allocate memory for particles
    Particle *particles = (Particle*)malloc(N * sizeof(Particle));
    if (particles == NULL) {
        printf("Error: Memory allocation failed.\n");
	free(particles);
        return 1;
    }
    
    // Read initial conditions from file
    FILE *input_file = fopen(filename, "rb");
    if (input_file == NULL) {
        printf("Error: Unable to open initial conditions file.\n");
	fclose(input_file);
        free(particles);
        return 1;
    }
    for (int i = 0; i < N ; i++) {
        fread(&particles[i].x, sizeof(double), 1, input_file);
        fread(&particles[i].y, sizeof(double), 1, input_file);
        fread(&particles[i].mass, sizeof(double), 1, input_file);
        fread(&particles[i].vx, sizeof(double), 1, input_file);
        fread(&particles[i].vy, sizeof(double), 1, input_file);
        fread(&particles[i].brightness, sizeof(double), 1, input_file);
        
        // Print values for debugging
        //printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, particles[i].brightness);
        //fprintf(input_file, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf", &particles[i].x, &particles[i].y; &particles[i].mass, &particles[i].vx, &particles[i].vy, &particles[i].brightness);
        
    }
    
    fclose(input_file);
    
    // Parallelize Simulation loop using OpenMP
    #pragma omp parallel for num_threads(num_threads) {
    for (int step = 0; step < nsteps; step++) {
       //split workload among threads
       #pragma omp for
       /*{
	  //Define thread specific arguments
	  ThreadArgs args;
          args.particles = particles;
          args.N = N;
	  args.start = 0;
	  args.end = N;
	  int tid = omp_get_thread_num();
	  int chunk_size = N / num_threads;
	  args.start = tid * chunk_size;
	  args.end = (tid == num_threads - 1) ? N : args.start + chunk_size;
	 */
       for(int i = 0; i < N; i++) {  
	  //Execute thread function
	  //simulation_force((void*)&i, particles, N, delta_t);
	  simulation_force((void*)&delta_t);
       }
      //Update particles position (outside thread function)
      #pragma omp barrier
      #pragma omp for
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
      return 1;
   }
   for (int i = 0; i < N; i++) {
        fwrite(&particles[i].x, sizeof(double), 1, output_file);
        fwrite(&particles[i].y, sizeof(double), 1, output_file);
        fwrite(&particles[i].mass, sizeof(double), 1, output_file);
        fwrite(&particles[i].vx, sizeof(double), 1, output_file);
        fwrite(&particles[i].vy, sizeof(double), 1, output_file);
        fwrite(&particles[i].brightness, sizeof(double), 1, output_file);
        
        // Print values for debugging
        //printf("Particle %d: x=%0.6lf, y=%0.6lf, mass=%0.6lf, vx=%0.6lf, vy=%0.6lf, brightness=%0.6lf\n", i, particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, particles[i].brightness);
        
   }
    
   fclose(output_file);
    
   // Free allocated memory
   free(particles);
   pthread_mutex_destroy(&mutex);
    return 0;
}

