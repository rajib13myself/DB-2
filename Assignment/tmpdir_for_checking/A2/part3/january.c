#include <stdio.h>
#include <stdlib.h>

//Designed a Stucture to preseve a day's temperature data.
struct DayTempData {
	int index;
	double minTemp;
	double maxTemp;
	struct DayTempData* next;
};

//Fuction to insert a day's Temp Data into the linked list
void insertDayTData(struct DayTempData** head, int index, double minTemp, double maxTemp) {
	struct DayTempData* current = *head;
	while (current != NULL) {
		if(current->index == index) {
			//Data for this day already exists, updating with new one.
			//printf("Data for day %d already exists. Ignoring dulciate and update with new one.\n", index);
			current->minTemp = minTemp;
			current->maxTemp = maxTemp;
			return;
		}
		current = current->next;
	}
	//If data doesn't exist, proceed with insertion.
	struct DayTempData* newNode = (struct DayTempData*)malloc(sizeof(struct DayTempData));
	if (newNode == NULL) {
		fprintf(stderr, "Memory Allocation error\n");
		exit(EXIT_FAILURE);
	}
	
	newNode->index = index;
	newNode->minTemp = minTemp;
	newNode->maxTemp = maxTemp;
	newNode->next = NULL;
	
	//Find the appropiate location to insert in a sorted  manner.
	//struct DayTempData* current = *head;
	current = *head;
	struct DayTempData* prev = NULL;
	
	//Data Insert with pointer movement into the linked list.
	while (current != NULL && current->index < index) {
		prev = current;
		current = current->next;
	}
	
	if (prev == NULL) {
		//Insert at the begaining position 
		newNode->next = *head;
		*head = newNode;
	} else {
		//Insert in the middle or at the end
		prev->next = newNode;
		newNode->next = current;
	}
}

//Fucntion to delete a day's data from the linked list
void deleteTempData(struct DayTempData** head, int index) {
	struct DayTempData* current = *head;
	struct DayTempData* prev = NULL;
	
	//Finding pointer 
	while (current != NULL && current->index != index) {
		prev =current;
		current = current->next;
	}
	if (current == NULL) {
	//	printf("Day with index %d not found\n", index);
		return;
	}
	
	if (prev == NULL) {
		//Delete from the beginning
		*head = current->next;
	} else {
		//Delete from the middle or end
		prev->next = current->next;
	}
	free(current);
}

//Function to point all data in the linked list
void printData(struct DayTempData* head) {
	printf("day      min         max\n");
	struct DayTempData* current = head;
	while (current != NULL) {
		printf("%d      %.6lf    %.6lf\n", current->index, current->minTemp, current->maxTemp);
		current = current->next;
	}
}

//Function to free the memory allocated for the linked list
void freeLinkedList(struct DayTempData* head) {
	struct DayTempData* current = head;
	while (current != NULL) {
		struct DayTempData* next = current->next;
		free(current);
		current = next;
	}
}

//Main function to run the Temperature program
int main() {
	struct DayTempData* database = NULL;
	char command;
	int index;
	double minTemp, maxTemp;
	
	while (1) {
		printf("Enter command: ");
		scanf(" %c", &command);
	
	switch (command) {
		case 'A':
			scanf("%d %lf %lf", &index, &minTemp, &maxTemp);
			insertDayTData(&database, index, minTemp, maxTemp);
			break;
		case 'D':
			scanf("%d", &index);
			deleteTempData(&database, index);
			break;
		case 'P':
			printData(database);
			break;
		case 'Q':
			freeLinkedList(database);
			return 0;
			break;
		default:
			printf("Invalid command\n");
			break;
		}
	}
	return 0;
}

	
