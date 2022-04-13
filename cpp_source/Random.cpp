#include "Random.h"

std::function<double()> default_double_generator(const double& min, const double& max)
{
	return [&min, &max]()
	{
		static std::random_device rd{};
		static std::mt19937	      gen{ rd() };
		static std::uniform_real_distribution<double> dist(min, max);

		return dist(gen);
	};
}

std::function<double()> default_int_generator(const double& min, const double& max)
{
	return [&min, &max]()
	{
		static std::random_device rd{};
		static std::mt19937	      gen{ rd() };
		static std::uniform_real_distribution<double> dist(min, max);

		return static_cast<int>(dist(gen));
	};
}
