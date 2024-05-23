#include <stdio.h>
#include <math.h>

int main() {
    // Test for becoming 'inf'
    float f = 1.0f;
   // printf("Float becoming 'inf':\n");
    while (isfinite(f)) {
      //  printf("%f\n", f);
        f *= 100.0f;
    }
   // printf("Float overflowed to 'inf'\n\n");

    // Test for producing 'nan'
    double d = sqrt(-1.0);
    printf("Double producing 'nan': %f\n\n", d);

    // Test operations on 'inf' and 'nan'
    double inf = INFINITY;
    double nan = NAN;
    printf("Operations on 'inf' and 'nan':\n");
    printf("inf + 1 = %f\n", inf + 1.0);
    printf("nan + 1 = %f\n\n", nan + 1.0);

    // Test relative precision of floating-point numbers
    double epsilon = 1.0;
    printf("Relative precision test:\n");
    while ((1.0 + epsilon) > 1.0) {
        epsilon *= 0.5;
    }
    printf("Smallest epsilon such that 1.0 + epsilon > 1.0: %e\n", epsilon);

    return 0;
}

