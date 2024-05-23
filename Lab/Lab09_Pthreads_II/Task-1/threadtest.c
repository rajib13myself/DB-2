#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

void* the_thread_func(void* arg) {
  /* Do something here? */
  
 int *a = malloc(sizeof(int)); //allocate memory for a variable
  if( a == NULL) {
   printf("Memory allocation failed.\n");
   pthread_exit(NULL);
  } 
  int prod = 7;
  printf("Enter a value for (a): ");
  scanf("%d", a);
  *a = *a * prod;
  return (void*)a;
}

int main() {
  printf("This is the main() function starting.\n");

  /* Start thread. */
  pthread_t thread;
  printf("the main() function now calling pthread_create().\n");
  if(pthread_create(&thread, NULL, the_thread_func, NULL) != 0) {
    printf("ERROR: pthread_create failed.\n");
    return -1;
  }

  printf("This is the main() function after pthread_create()\n");

  /* Do something here? */
  void* output;
  /* Wait for thread to finish. */
  printf("the main() function now calling pthread_join().\n");
  if(pthread_join(thread, &output) != 0) {
    printf("ERROR: pthread_join failed.\n");
    return -1;
  }
  printf("output of function thread: %d\n", *(int*)output);
  free(output); //Free the thread fucntion memory.

  return 0;
}
