#pragma once

#include <chrono>

class Timer
{


public:

	void start ();
	void stop  ();

	double elapsed();


private:

	long long														_duration = 0;
	std::chrono::time_point<std::chrono::high_resolution_clock>		_start_point;

};