#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *input_file = fopen("result.gal", "rb"); // Open binary input file
    if (input_file == NULL) {
        perror("Error opening input file");
        return 1;
    }

    FILE *output_file = fopen("output_re.txt", "w"); // Open text output file
    if (output_file == NULL) {
        perror("Error opening output file");
        fclose(input_file);
        return 1;
    }

    // Read and write data until end of input file is reached
    int value;
    while (fread(&value, sizeof(int), 1, input_file) == 1) {
        fprintf(output_file, "%d\n", value); // Write integer value to text file
    }

    // Close both files
    fclose(input_file);
    fclose(output_file);

    printf("Conversion successful!\n");

    return 0;
}

