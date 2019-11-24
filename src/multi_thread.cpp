#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "multi_thread.h"
#include <pthread.h>


void *wait(void *t) {
   long tid = (long)t;
   pthread_exit(NULL);
}

int run_task(int total, int num_cores, void* task(void* data), void* data) {
   int ret;
   int i;
   pthread_t* threads = (pthread_t*) malloc(total* sizeof(pthread_t));
   pthread_attr_t attr;
   void *status;
	 int task_num = 0;
	 for( i = 0; i < total; i+= num_cores) {
		  task_num = std::min(num_cores, total-num_cores);
			for (int j = 0; j < task_num; ++j) {
				thread_data* thd;
				thd->thread_id = i+j;
			  thd->data  = data;
				ret = pthread_create(&threads[i+j], &attr, task, thd);
				if (ret) {
					 printf("Error:unable to create thread %d err code %d\n", i+j, ret);
					 return -1;
				}
			}
		 for(int j = 0; j < task_num; j++) {
				ret = pthread_join(threads[i+j], &status);
				if (ret) {
					 printf("Error:unable to join thread %d\n", i);
					 return -1;
				}
		 }
   }

	free(threads);
	return 0;
}
