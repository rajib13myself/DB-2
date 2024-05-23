#include <stdio.h>
#include <stdlib.h>

typedef struct product {
    char name[50];
    double price;
} product_t;

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;  // Exit with an error code
    }

    FILE *file;
    product_t *arr_of_prod;
    int totalProducts;

    // Open the file for reading
    file = fopen(argv[1], "r");

    // Check if the file is successfully opened
    if (file == NULL) {
        printf("Could not open the file.\n");
        return 1;  // Exit with an error code
    }

    // Read the total number of products from the first line
    if (fscanf(file, "%d", &totalProducts) != 1) {
        printf("Error reading the total number of products.\n");
        fclose(file);
        return 1;  // Exit with an error code
    }

    // Dynamically allocate memory for the array of products
    arr_of_prod = (product_t *)malloc(totalProducts * sizeof(product_t));

    // Check if memory allocation is successful
    if (arr_of_prod == NULL) {
        printf("Memory allocation failed.\n");
        fclose(file);
        return 1;  // Exit with an error code
    }

    // Read and store each product in the dynamically allocated array
    for (int i = 0; i < totalProducts; i++) {
        if (fscanf(file, "%s %lf", arr_of_prod[i].name, &arr_of_prod[i].price) != 2) {
            printf("Error reading product data.\n");
            free(arr_of_prod);
            fclose(file);
            return 1;  // Exit with an error code
        }
    }

    // Output the table header
    printf("%-10s %-10s\n", "Product", "Price");
    printf("--------------------\n");

    // Output each product from the dynamically allocated array
    for (int i = 0; i < totalProducts; i++) {
        printf("%-10s %-10.2lf\n", arr_of_prod[i].name, arr_of_prod[i].price);
    }

    // Free the dynamically allocated memory
    free(arr_of_prod);

    // Close the file
    fclose(file);

    return 0;
}

