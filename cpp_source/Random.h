#pragma once
#include <random>
#include <functional>

std::function<double()> default_double_generator (const double& min, const double& max);
std::function<double()> default_int_generator	 (const double& min, const double& max);