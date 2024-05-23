#include <stdio.h>
#include <math.h>

int main() {
	int choice;
	//Prompt user choice.
	printf("Choose an option:\n");
	printf("1: Sum or Product of Two Integers.\n");
	printf("2: Find Largest value from three real numbers.\n");
	printf("3: Find Second Largest from three real numbers.\n");
	scanf("%d", &choice);
	
	switch (choice) {
		case 1: {
			// Sum or Product of Two Integers.
			//declare variable for 2 integer numbers.
			int num1, num2;
			printf("\nPlease Enter two integer numbers: ");
			scanf("%d %d", &num1, &num2);
			//Check two variables even or not, if both are even then sum otherwise product.
			//int sum, prd;
			if((num1 % 2 == 0) && (num2 % 2 == 0)) {
				printf("\nBoth input numbers are even and sum of the numbers is : %d\n", num1 + num2);
			} else {
				printf("\nBoth input numbers are not even and product of the numbers is : %d\n", num1 * num2);
			}
			break;
		}
		case 2: {
			//Largest Number by Absolute Value.
			double num1, num2, num3;
			//Read Three real numbers
			printf("\nEnter Three real numbers betwen space of them.: ");
			scanf("%lf %lf %lf", &num1, &num2, &num3);
			//Find the largest from input values.
			double abs1 = fabs(num1);
			double abs2 = fabs(num2);
			double abs3 = fabs(num3);
			
			double maxabs = fmax(abs1, fmax(abs2, abs3));
			//Output of the largest value among three real numbers.
			printf("\nThe Largest absolute value : %.2lf\n", maxabs);
			break;
		}
		case 3: {
			//Sceond Largest absolute value.
			double num1, num2, num3;
			//Read Three real numbers
			printf("\nEnter Three real numbers between space of them.: ");
			scanf("%lf %lf %lf", &num1, &num2, &num3);
			//Find the 2nd Largest among three real numbers.
			double abs1 = fabs(num1);
			double abs2 = fabs(num2);
			double abs3 = fabs(num3);
			double maxabs = fmax(abs1, fmax(abs2, abs3));
			double secondmax;
			if(maxabs == abs1) {
				secondmax = fmax(abs2, abs3);
			} else if (maxabs == abs2) {
				secondmax = fmax(abs1, abs3);
			} else {
				secondmax = fmax(abs1, abs2);
			}
			printf("\nThe Second Largest Absolute value is : %.2lf\n", secondmax);
			break;
		}
		default:
			printf("Invalid choice\n");
			break;
	}
	return 0;
}

