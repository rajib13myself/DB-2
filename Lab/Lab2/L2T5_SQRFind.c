#include <stdio.h>
#include <math.h>

int main() {

	//Find the perfect square root of input value.
	//Declare variable.
	int num1;
	//read value for variable.
	printf("\nEnter value for perfect square root finding.\n");
	scanf("%d", &num1);
	//Indentify input value is Square root or not.
	if(num1 >0) {
		double squareroot = sqrt(num1);
		int perfectsqrt = (int)squareroot;
		if((squareroot - perfectsqrt) == 0) {
			printf("%d is perfect square.\n", num1);
		} else {
			printf("%d is not perfect square.\n", num1);
		}
	}
	return 0;
}
