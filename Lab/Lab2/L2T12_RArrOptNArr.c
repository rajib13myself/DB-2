#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool is_prime(int num) {
	if(num <=1) {
		return false;
	}
	for(int i = 2; i * i <= num; i++) {
		if(num % i == 0) {
			return false;
		} 
	}
	return true;
}


int main() {
	//int *Originalarr = NULL;	//Pointer of an array
	int n = 0;		//Size of the array
	//int input = 0;		//Variable for value input
	//int *Newarr = NULL;	//Pointer of new Array
	//int size = 0;
	printf("Enter a value for array size: ");
	scanf("%d", &n);
	
	int *Originalarr = (int *)malloc(n * sizeof(int)); //allocate memory
	int OriginalSize = 0;
	//Read n integer value in the array.
	printf("\nEnter %d integer numbers: ", n);
	for(int i = 0; i < n; i++) {
		scanf("%d", &Originalarr[i]);
		OriginalSize++;
	}
	//Allocate new array
	int *newArray = NULL;
	int newSize = 0;
	for(int i = 0; i < OriginalSize; i++) {
		 if(is_prime(Originalarr[i])) {
			//Reallocate memory for new array.
			newArray = (int *)realloc(newArray, (newSize + 1) * sizeof(int));
			//Check Memory allocation is successfull
			if(newArray == NULL) {
				printf("Memory Allocation failed.\n");
				free(Originalarr);
				return 1; //Exit with an error code.
			}
			//Store prime number in new array.
			newArray[newSize] = Originalarr[i];
			newSize++;
		}
	}
		//print the element of new array.
		printf("Elements of the new Array: ");
		for(int i = 0; i <newSize; i++) {
			printf("%d ", newArray[i]);
		} 
		printf("\n");
		//Print size of new array,
		printf("Size of the new Array : %d \n", newSize);
		//Free the allocated memory.
		free(Originalarr);
		free(newArray);
		
				
	return 0;
}
