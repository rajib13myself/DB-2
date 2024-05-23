#include <stdio.h>
#include <stdlib.h>

// Function to calculate the factorial of a number
long long factorial(int n) {
    if (n == 0 || n == 1) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}

// Function to calculate the binomial coefficient C(n, k)
long long binomialCoefficient(int n, int k) {
    return factorial(n) / (factorial(k) * factorial(n - k));
}

// Function to print Pascal's Triangle
void printPascalsTriangle(int numRows) {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j <= i; j++) {
            printf("%lld ", binomialCoefficient(i, j));
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {

        // Check if the correct number of command-line arguments is provided
        if (argc != 2) {
                printf("Usage: %s <number_of_rows>\n", argv[0]);
                return 1;
        }
        // Get the number of rows from the command-line argument
        int numRows = atoi(argv[1]);
        // Check if the input is valid
        if (numRows < 0) {
                printf("Number of rows should be non-Negetive.\n");
                return  1;
        }
        // Print Pascal's Triangle
        printPascalsTriangle(numRows);

        return 0;
}
