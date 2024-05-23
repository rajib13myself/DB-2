#include <stdio.h>
#include <stdlib.h>

// Initial size of the array chunk
#define CHUNK_SIZE 10

int main() {
    int *arr = NULL;  // Pointer to the array
    int size = 0;     // Size of the array
    int capacity = CHUNK_SIZE; // Initial capacity of the array
    int input;
    int sum = 0;

    arr = (int *)malloc(capacity * sizeof(int));  // Allocate initial memory

    // Check if memory allocation is successful
    if (arr == NULL) {
        printf("Memory allocation failed.\n");
        return 1;  // Exit with an error code
    }
    printf("Input: ");
    while (scanf("%d", &input) == 1) {
        // Store the input in the array
        arr[size] = input;

        // Increment the size of the array
        size++;

        // Add the input to the sum
        sum += input;

        // Check if the array needs to be resized
        if (size == capacity) {
            // Increase the capacity by CHUNK_SIZE
            capacity += CHUNK_SIZE;

            // Reallocate memory for the array
            int *temp = (int *)realloc(arr, capacity * sizeof(int));

            // Check if memory allocation is successful
            if (temp == NULL) {
                printf("Memory allocation failed.\n");
                free(arr);
                return 1;  // Exit with an error code
            }

            arr = temp;
        }
    }

    // Print entered numbers
    printf("Entered numbers: ");
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    // Print the sum of entered numbers
    printf("Sum: %d\n", sum);

    // Free the allocated memory
    free(arr);

    return 0;
}
 	


	
