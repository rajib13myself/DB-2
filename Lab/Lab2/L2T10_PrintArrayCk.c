#include <stdio.h>
#include <stdlib.h>  // Include for malloc and free

// Function to print all values of an integer array
void print_array(int arr[], int size) {
    printf("Array values: ");
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

int main() {
	int *arr;
    	int n;
	printf("Enter a value for n: ");
    scanf("%d\n", &n);

    // Allocate memory for the array
    arr = (int *)malloc(n * sizeof(int));

    // Check if memory allocation is successful
    if (arr == NULL) {
        printf("Memory allocation failed.\n");
        return 1;  // Exit with an error code
    }

    // Fill the array with random numbers
    for (int i = 0; i < n; ++i) arr[i] = rand() % 100;  // random number from 0 to 99

    // Call the function to print the array
    print_array(arr, n);

    // Free the allocated memory
    free(arr);

    return 0;
}

