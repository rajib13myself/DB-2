#include <stdio.h>
#include <math.h>

// Define a structure for a 2D vector
typedef struct {
    double x;
    double y;
} Vector2D;

// Function to subtract two vectors
Vector2D subtract(Vector2D v1, Vector2D v2) {
    Vector2D result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    return result;
}

// Function to calculate the magnitude of a vector
double magnitude(Vector2D v) {
    return sqrt(v.x * v.x + v.y * v.y);
}

// Function to normalize a vector
Vector2D normalize(Vector2D v) {
    double mag = magnitude(v);
    Vector2D result;
    result.x = v.x / mag;
    result.y = v.y / mag;
    return result;
}

int main() {
    // Define two vectors xi, xj, yi, and yj
    Vector2D xi = {1.0, 2.0};  // Example values
    Vector2D xj = {3.0, 4.0};  // Example values
    Vector2D yi = {5.0, 6.0};  // Example values
    Vector2D yj = {7.0, 8.0};  // Example values

    // Calculate rij
    Vector2D rij = subtract(xi, xj);

    // Calculate rij squared
    double r2_ij = pow(rij.x, 2) + pow(rij.y, 2);

    // Calculate normalized rij
    Vector2D normalized_rij = normalize(rij);

    // Output results
    printf("rij = (%.2f, %.2f)\n", rij.x, rij.y);
    printf("r2_ij = %.2f\n", r2_ij);
    printf("^rij = (%.2f, %.2f)\n", normalized_rij.x, normalized_rij.y);

    return 0;
}

