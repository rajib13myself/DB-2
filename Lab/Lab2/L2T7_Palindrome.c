#include <stdio.h>
#include <math.h>

int main() {
	//Palindrome Check.
	//declare variables.
	int number, originalNum, reverseNum = 0, reminder;
	printf("Enter a number: ");
	scanf("%d", &number);
	//Save original number.
	originalNum = number;
	//Verify Input number is Palindrome  or not.
	while (number !=0) {
		reminder = number % 10;
		reverseNum = reverseNum * 10 + reminder;
		number /=10;
	}
	if(originalNum == reverseNum) {
		printf("\nIt is a Palindrome.\n");
	} else {
		printf("\nIt is not a Palindrome.\n");
	}
	return 0;
}
