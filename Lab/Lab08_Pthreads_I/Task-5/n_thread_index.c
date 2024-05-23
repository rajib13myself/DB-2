#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

// Structure for thread data
typedef struct {
    int thread_index;
} ThreadData;

// Thread function
void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    printf("Thread %d\n", data->thread_index);
    pthread_exit(NULL);
}

int main() {
    int num_threads;
    printf("Enter the number of threads (N): ");
    scanf("%d", &num_threads);

    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];

    // Create threads
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].thread_index = i;
        pthread_create(&threads[i], NULL, thread_function, (void*)&thread_data[i]);
    }

    // Join threads
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    return 0;
}

