#include "Timer.h"

void Timer::start()
{
	_start_point = std::chrono::high_resolution_clock::now();
}

void Timer::stop()
{
	auto end_point = std::chrono::high_resolution_clock::now();

	auto start = std::chrono::time_point_cast<std::chrono::microseconds>(_start_point).time_since_epoch().count();
	auto end   = std::chrono::time_point_cast<std::chrono::microseconds>(end_point).time_since_epoch().count();

	_duration = end - start;
}

double Timer::elapsed()
{
	return _duration * 0.001;
}
