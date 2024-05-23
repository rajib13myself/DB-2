#include <stdio.h>

int main() {
    FILE *file;
    char product[50];
    double price;

    // Open the file for reading
    file = fopen("/home/rdatta3822/HPP/Lab/Lab2/data.txt", "r");

    // Check if the file is successfully opened
    if (file == NULL) {
    perror("Error opening file");
    return 1;  // Exit with an error code
	}


    // Read data from the file and output as a table
    printf("%s %s\n", "Product", "Price");
    printf("--------------------\n");

    while (fscanf(file, "%s %lf", product, &price) == 2) {
        printf("%s %.2lf\n", product, price);
    }

    // Close the file
    fclose(file);

    return 0;
}

