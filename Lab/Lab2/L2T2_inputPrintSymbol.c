#include <stdio.h>
#include <ctype.h>

int main() {
	printf("\nPrint Symbol with orientation by user input.");
	//declare integer variable
	int a, b;
	printf("\nEnter a number value for a variable: ");
	scanf("%d", &a);
	printf("\nEnter a number value for b variable: ");
	scanf("%d", &b);
	//Output the rectangle using '*' and '.' symbols
			
	//if (isdigit((char)a) && isdigit((char)b)) {
		for(int i = 1; i <= a; i++) {
			for(int j =1; j <= b; j++) {
				//Ouput '*'  for the border and '.' for the interior 
				if (i == 1 || i == a || j == 1 || j == b) {
					printf("*");
				} else {
					printf(".");
				}
			}
			//Move to the next line after completing a row.
			printf("\n");
		}
		return 0;
	//}

}
