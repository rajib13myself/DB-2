#include <stdio.h>
#include <time.h> // Include time.h for clock_t and clock

#define ITERATIONS 100000000

int add(int a, int b) {
    return a + b;
}

inline int add_inline(int a, int b) {
    return a + b;
}

int main() {
    clock_t start_non_inline, end_non_inline, start_inline, end_inline;
    double time_non_inline, time_inline;
    int result = 0;

    // Measure time for non-inlined function
    start_non_inline = clock();
    for (int i = 0; i < ITERATIONS; ++i) {
        result = add(i, i);
    }
    end_non_inline = clock();
    time_non_inline = ((double)(end_non_inline - start_non_inline)) / CLOCKS_PER_SEC;

    // Measure time for inlined function
    start_inline = clock();
    for (int i = 0; i < ITERATIONS; ++i) {
        result = add_inline(i, i);
    }
    end_inline = clock();
    time_inline = ((double)(end_inline - start_inline)) / CLOCKS_PER_SEC;

    // Display results
    printf("Time taken for non-inlined function: %f seconds\n", time_non_inline);
    printf("Time taken for inlined function: %f seconds\n", time_inline);

    return 0;
}

