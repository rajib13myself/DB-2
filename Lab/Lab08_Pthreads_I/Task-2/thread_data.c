#include <stdio.h>
#include <pthread.h>

void* the_thread_func(void* arg) {
  /* Do something here? */
  printf("the_thread_func is now starting to fucntion here.\n");
  int i;
  double multiply = 1.0;
  for(i = 1; i<9999; i++)
    multiply *= 1.000007;
  printf("Multipy value for the thread function: %0.3f\n", multiply);
  return NULL;
}

int main() {
  printf("This is the main() function starting.\n");

  double data_for_thread[3];
  data_for_thread[0] = 5.7;
  data_for_thread[1] = 9.2;
  data_for_thread[2] = 1.6;
  double data_for_thread_B[3];
  data_for_thread_B[0] = 1.7;
  data_for_thread_B[1] = 3.7;
  data_for_thread_B[2] = 5.7;




  /* Start thread. */
  pthread_t thread;
  printf("the main() function now calling pthread_create().\n");
  pthread_create(&thread, NULL, the_thread_func, data_for_thread);

  printf("This is the main() function after pthread_create()\n");

  /* Do something here? */
  int i;
  double prod = 1.0;
  double total_prod = 1.0;
  char strSentence[27];
  char current_letter = 'a';
  for(i = 0; i<26; i++) {
    strSentence[i] = current_letter;
    current_letter++;
    total_prod = prod * i;
  }
  // Null-terminate the string
  strSentence[26] = '\0';

  // Print the resulting string
  printf("Letters from 'a' to 'z': %s\n", strSentence);
  //Print production stering;
  printf("Product of the total_production: %0.3f\n", total_prod);
 // return 0;
  /* Wait for thread to finish. */
  printf("the main() function now calling pthread_join().\n");
  pthread_join(thread, NULL);

  return 0;
}
