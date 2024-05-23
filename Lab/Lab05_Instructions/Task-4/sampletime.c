#include <stdio.h>
#include <time.h>

#define ITERATIONS 1000000000 // Adjust the number of iterations as needed

int main() {
    unsigned long i;
    double result = 1.0; // To prevent compiler optimization

    // Measure time for addition
    clock_t start_addition = clock();
    for (i = 0; i < ITERATIONS; ++i) {
        result += i;
    }
    clock_t end_addition = clock();
    double time_addition = ((double) (end_addition - start_addition)) / CLOCKS_PER_SEC;

    // Measure time for subtraction
    clock_t start_subtraction = clock();
    for (i = 0; i < ITERATIONS; ++i) {
        result -= i;
    }
    clock_t end_subtraction = clock();
    double time_subtraction = ((double) (end_subtraction - start_subtraction)) / CLOCKS_PER_SEC;

    // Measure time for multiplication
    clock_t start_multiplication = clock();
    for (i = 1; i <= ITERATIONS; ++i) {
        result *= i;
    }
    clock_t end_multiplication = clock();
    double time_multiplication = ((double) (end_multiplication - start_multiplication)) / CLOCKS_PER_SEC;

    // Measure time for division
    clock_t start_division = clock();
    for (i = 1; i <= ITERATIONS; ++i) {
        result /= i;
    }
    clock_t end_division = clock();
    double time_division = ((double) (end_division - start_division)) / CLOCKS_PER_SEC;

    // Display results
    printf("Time taken for addition: %f seconds\n", time_addition);
    printf("Time taken for subtraction: %f seconds\n", time_subtraction);
    printf("Time taken for multiplication: %f seconds\n", time_multiplication);
    printf("Time taken for division: %f seconds\n", time_division);

    return 0;
}

