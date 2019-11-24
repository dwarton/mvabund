#include <ctime>
#include <iostream>
#include <string>
#include <stdlib.h>

#ifndef  TIME_PROFILE_H_
#define TIME_PROFILE_H_
class TimeDebug {
public:
  TimeDebug() : start_(std::clock()) {}
  void start() { start_ = std::clock(); }
  float elapsed_time() {
  	end_ = std::clock();
		//printf("elapsed_time %f seconds\n", (float(end_ - start_))/CLOCKS_PER_SEC);
	  return static_cast<float>(end_ - start_) / CLOCKS_PER_SEC;
	}
	std::string currentTime() {
			time_t now = time(NULL); 
			tm *ltm = localtime(&now); 
			std::string var;
			var += std::to_string(1900 + ltm->tm_year);
			var += "-" + std::to_string(1 + ltm->tm_mon);
   	  var += "-" + std::to_string(ltm->tm_mday);
			var += " " + std::to_string(ltm->tm_hour);
			var += ":" + std::to_string(ltm->tm_min);
			var += ":" + std::to_string(ltm->tm_sec);
		  //printf("\nCurrent time is %s\n", var);
			return var;
	}

private:
  clock_t start_;
  clock_t end_;
};
#endif
