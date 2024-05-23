#include <stdio.h>
#include <math.h>

int main() {
	//Compute Quotient and Reminder for input values division.
	//declare variables.
	int divident, divisor;
	//Read values for variables;
	printf("\nEnter divident: ");
	scanf("%d", &divident);
	printf("\nEnter divisor: ");
	scanf("%d", &divisor);
	//Indentify division fucntion and its result.
	if((divident % divisor) == 0) {
		int quotient = divident / divisor;
		int reminder = (divident % divisor);
		printf("\nQuotient: %d\n", quotient);
		printf("Reminder: %d\n", reminder);
	} else {
		int quotient = (divident / divisor);
		int reminder = (divident % divisor);
		printf("\nQuotient: %d\n", quotient);
		printf("\nReminder: %d\n", reminder);
	}
	return 0;
}

