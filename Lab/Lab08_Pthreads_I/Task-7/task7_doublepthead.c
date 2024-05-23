#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

#define NUM_THREADS 2

void* thread_function(void* arg) {
    int thread_index = *(int*)arg;
    printf("Thread %d created by Thread %lu\n", thread_index, pthread_self());

    // Create two new threads
    pthread_t sub_threads[2];
    for (int i = 0; i < 2; ++i) {
        int new_thread_index = thread_index * 2 + i + 1; // Calculate new thread index
        pthread_create(&sub_threads[i], NULL, thread_function, &new_thread_index);
    }

    // Wait for the sub-threads to finish
    for (int i = 0; i < 2; ++i) {
        pthread_join(sub_threads[i], NULL);
    }

    printf("Thread %d finished\n", thread_index);
    pthread_exit(NULL);
}

int main() {
    pthread_t threads[NUM_THREADS];

    // Create two threads
    for (int i = 0; i < NUM_THREADS; ++i) {
        int thread_index = i + 1;
        pthread_create(&threads[i], NULL, thread_function, &thread_index);
    }

    // Wait for the threads to finish
    for (int i = 0; i < NUM_THREADS; ++i) {
        pthread_join(threads[i], NULL);
    }

    printf("Main thread finished\n");
    return 0;
}

