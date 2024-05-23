#include <stdio.h>
#include <pthread.h>

void* the_thread_func(void* arg) {
  /* Do something here? */
  printf("the_thread_func() starting doing something here\n");
  long int i;
  double sub;
  for(i = 0; i<100000000; i++)
      sub -= 1;
  printf("The result of work in the thread function(): sub=%f\n", sub);
  return NULL;
}

int main() {
  printf("This is the main() function starting.\n");

  /* Start thread. */
  pthread_t thread;
  printf("the main() function now calling pthread_create().\n");
  pthread_create(&thread, NULL, the_thread_func, NULL);

  printf("This is the main() function after pthread_create()\n");

  /* Do something here? */
  printf("main() starting doing some work.\n");
  long int i;
  double sum;
  for(i = 0; i < 1000000000; i++){
   sum += 1e-7;
} 

  printf("Result of work in main(): sum = %f\n", sum); 
   printf("Print loop number : %ld\n", i);

  /* Wait for thread to finish. */
  printf("the main() function now calling pthread_join().\n");
  pthread_join(thread, NULL);

  return 0;
}
