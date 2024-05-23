#include <stdio.h>

int main(void) {
	FILE *file = fopen("little_bin_file", "rb"); //open the file in binary mode.
	if (file == NULL) {
		perror("Error opening file"); 	//Error message if file has null character.
		return 1;
	}
	
	//declare four types of variable
	int integerData;	
	double doubleData;
	char charData;
	float floatData;

	//Read Data from the file.
	fread(&integerData, sizeof(int), 1, file);
	fread(&doubleData, sizeof(double), 1, file);
	fread(&charData, sizeof(char), 1, file);
	fread(&floatData, sizeof(float), 1, file);
	
	//Close the file
	fclose(file);	

	//size_t total_DataSize;	//declare variable for total data size
	
	//Calulate the sum of binary data size
	//total_DataSize = sizeof(integerData) + sizeof(doubleData) + sizeof(charData) + sizeof(floatData);

	//Print the Data  and the sum of binary Data size.
	if(integerData == (int)integerData) {
		//printf("\nAn integer number represented using the datatype int : %d", integerData);
		printf("%d\n", integerData);
	} if(doubleData == (double)doubleData) {
		//printf("\nA foating-point number represented using the datatype double: %lf", doubleData);
		printf("%lf\n", doubleData);
	} if(charData == (char)charData) {
		//printf("\nA character represented using the datatype char: %c", charData);
		printf("%c\n", charData);
	} if(floatData == (float)floatData) {
		//printf("\nA foating-point number represented using the datatype float: %.1lf\n", floatData);
		printf("%.1lf\n", floatData);
	}
	//printf("\nSum of the binary data sizes: \033[0m %zu bytes\n", total_DataSize);
	
	return 0;
	
}

