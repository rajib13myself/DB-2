#include <stdio.h>
#include <math.h>

#define G 6.674e-11 // gravitational constant
#define EPSILON 1e-3 // small number for stability
#define DT 0.01 // time step size

typedef struct {
    double x, y; // position
    double vx, vy; // velocity
    double ax, ay; // acceleration
    double mass; // mass
} Particle;


// Function to calculate distance between two particles
double distance(Particle p1, Particle p2) {
    return sqrt(((p1.x - p2.x) * (p1.x - p2.x))  + ((p1.y - p2.y) * (p1.y - p2.y)));
}

//Function for distance vector(Particle p1, Particle p2)
double dist_vector(Particle p1, Particle p2) {
    return ((p1.x - p2.x) * 1.0 + (p1.y - p2.y) * 1.0); //Here ex and ey are unit vector
}

//Function for normalize both particles distance
double norm_dist(Particle p1, Particle p2) {
    double r_vec = dist_vector(p1, p2);
    double r_dist = distance(p1, p2);
    double r_hat = r_vec / r_dist;
    return r_hat;
}
// Function to calculate force between two particles using modified force formula
double force(Particle p1, Particle p2) {
    double r_dist = distance(p1, p2);
    double r_vec = dist_vector(p1, p2);
    double f = -G * p1.mass * p2.mass * ( r_vec / ((r_dist + EPSILON) * (r_dist + EPSILON) * (r_dist + EPSILON)));
    return f;
}

int main() {
    // Number of particles
    int N = 3;
    
    // Array of particles
    Particle particles[N];
    
    // Initialize particles
    particles[0].x = 0;
    particles[0].y = 0;
    particles[0].vx = 0;
    particles[0].vy = 0;
    particles[0].mass = 1.0;
    
    particles[1].x = 1;
    particles[1].y = 0;
    particles[1].vx = 0;
    particles[1].vy = 0;
    particles[1].mass = 1.0;
    
    particles[2].x = 0;
    particles[2].y = 1;
    particles[2].vx = 0;
    particles[2].vy = 0;
    particles[2].mass = 1.0;
    
    // Time loop
    for (double t = 0; t < 10; t += DT) {
        // Calculate forces and update accelerations
        double f;
        for (int i = 0; i < N; i++) {
            particles[i].ax = 0;
            particles[i].ay = 0;
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    f = f + force(particles[i], particles[j]);
                    //double r = distance(particles[i], particles[j]);
                    //particles[i].ax += f * (particles[j].x - particles[i].x) / r;
                    //particles[i].ay += f * (particles[j].y - particles[i].y) / r;
                
                }
            }
            particles[i].ax += f * particles[i].mass;
            particles[i].ay += f * particles[i].mass; 
        }
        
        // Update velocities and positions
        for (int i = 0; i < N; i++) {
            if (i < N-1) {
                particles[i+1].vx = particles[i].vx + particles[i].ax * DT;    // / particles[i].mass;
                particles[i+1].vy = particles[i].vy + particles[i].ay * DT;    // particles[i].mass;
                particles[i+1].x = particles[i].x + particles[i+1].vx * DT;
                particles[i+1].y = particles[i].y + particles[i+1].vy * DT;
            }
            
        }
        
        // Output positions
        for (int i = 0; i < N; i++) {
            printf("Particle %d: x = %lf, y = %lf, mass = %lf\n", i+1, particles[i].x, particles[i].y, particles[i].mass);
        }
        printf("\n");
    }
    
    return 0;
}
