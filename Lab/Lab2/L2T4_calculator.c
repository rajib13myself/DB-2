#include <stdio.h>
#include <math.h>

int main() {
	//Simple Calculator Program
	printf("\nSimple Calculator Program.\n");
	//declare variables for input numbers and symbol for operation.
	int a, b;
	char symbol;
	//Read values for variables
	printf("\nPlease enter values and symbol('+'||'-'||'*'||'/'||'^') for calculation.\n");
	scanf("%d %c %d", &a, &symbol, &b);
	//char choice;
	//Indentify symbol for operation and calculation.
	switch (symbol) {
		case '+': {
			printf("Sum: %d\n", a + b);
			break;
		}
		case '-': {
			printf("Subtruction: %d\n", a - b);
			break;
		}
		case '*': {
			printf("Multiplication: %d\n", a * b);
			break;
		}
		case '/': {
			if(b != 0) {
				printf("Division: %d\n", a / b);
			} else {
				printf("Division result infinity as 2nd number is zero.\n");
			}
			break;
		}
		case '^': {
			printf("b times a is: %.2lf\n", pow(a, b));
			break;
		}	
		default:
			printf("Invalid symbol\n");
			break;
	}

}

