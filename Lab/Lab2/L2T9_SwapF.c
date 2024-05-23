#include <stdio.h>


//Function of Swap values of Two integers.
void swapNums(int *a, int *b) {
	int temp = *a;
	*a = *b;
	*b = temp;
}

//Function of Swap values of two pointers (string)
void swapPointers(char **s1, char **s2) {
	char *temp = *s1;
	*s1 = *s2;
	*s2 = temp;
}

int main() {
    // Declare and initialize two integers
    int a = 3, b = 4;

    // Declare and initialize two strings (pointers to char arrays)
    char *s1 = "first";
    char *s2 = "second";

    // Output original values
    printf("a=%d, b=%d\n", a, b);
    printf("s1=%s, s2=%s\n", s1, s2);

    // Swap integers using swapNums function
    swapNums(&a, &b);

    // Swap pointers using swapPointers function
    swapPointers(&s1, &s2);

    // Output values after swapping
    printf("a=%d, b=%d\n", a, b);
    printf("s1=%s, s2=%s\n", s1, s2);

    return 0;
}
