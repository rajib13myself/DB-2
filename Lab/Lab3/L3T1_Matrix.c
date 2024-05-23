#include <stdio.h>

void create2DMMatrix(int n) {
	//Need to verify n should be greater than 1.
	if(n <= 1) {
		printf("Please enter a value of n greater than 1. \n");
		return;
	}
	int matrixVal[n][n]; //Declare 2D Array
	
	//Intialize the matrix as per given pattern.
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < n; ++j) {
			if( i == j) {
				matrixVal[i][j] = 0;
			} else {
				matrixVal[i][j] = (i > j) ? -1 : 1;
			}
		}
	}
	//Display the Matrix.
	printf("Generated Matrix:\n");
	for(int i = 0; i < n; ++i) {
		printf("| ");
		for(int j = 0; j < n; ++j) {
			printf("%d  ", matrixVal[i][j]);
		}
		printf(" |\n");
	} 


}


int main () {
	//Declare a variable for matrix dimention n x n
	int n;
	printf("Enter a value for matrix size: ");
	scanf("%d", &n); 	//Read the value for matrix size.
	
	//Create and Display 2D Array as per given pattern.
	create2DMMatrix(n);
	
	return 0;
}
