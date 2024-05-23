#include <stdio.h>

int main() {
    FILE *file;
    char product[50];
    double price;
    int totalProducts;

    // Open the file for reading
    file = fopen("data.txt", "r");

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

    // Output the table header
    printf("%-10s %-10s\n", "Product", "Price");
    printf("--------------------\n");

    // Read and output each product and price
    for (int i = 0; i < totalProducts; i++) {
        if (fscanf(file, "%s %lf", product, &price) != 2) {
            printf("Error reading product data.\n");
            fclose(file);
            return 1;  // Exit with an error code
        }

        // Output the product and price
        printf("%-10s %-10.2lf\n", product, price);
    }

    // Close the file
    fclose(file);

    return 0;
}

