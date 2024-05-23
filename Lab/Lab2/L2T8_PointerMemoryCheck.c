#include <stdio.h>

int main() {
    // Declare and initialize variables
    double myDouble = 3.14;
    int myInt = 42;
    char myChar = 'A';

    // Output information for double
    printf("Double Value: %lf\n", myDouble);
    printf("Address: %p\n", (void*)&myDouble);
    printf("Memory Size: %lu bytes\n\n", sizeof(myDouble));

    // Output information for int
    printf("Integer Value: %d\n", myInt);
    printf("Address: %p\n", (void*)&myInt);
    printf("Memory Size: %lu bytes\n\n", sizeof(myInt));

    // Output information for char
    printf("Char Value: %c\n", myChar);
    printf("Address: %p\n", (void*)&myChar);
    printf("Memory Size: %lu bytes\n\n", sizeof(myChar));

    return 0;
}

