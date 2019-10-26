#include <iostream>
#include "random.h"


RandomNumbers::RandomNumbers(unsigned long int s)
:seed(s)
{
	if(seed == 0)
	{
		std::random_device rd;
		seed = rd();
	}
	rng = std::mt19937(seed);
}

void RandomNumbers::uniform_double(std::vector<double>& random_list, double lower, double upper)
{
	std::uniform_real_distribution<> d(lower,upper);
	for (size_t i(0); i<random_list.size(); ++i)
	{
		random_list[i] = d(rng);
	}
}

double RandomNumbers::uniform_double(double lower, double upper)
{
	std::uniform_real_distribution<> d(lower,upper);
	return d(rng);
}

void RandomNumbers::normal(std::vector<double>& random_list, double mean, double sd)
{
	std::normal_distribution<> d(mean,sd);
	for (size_t i(0); i<random_list.size(); ++i)
	{
		random_list[i] = d(rng);
	}
}

double RandomNumbers::normal(double mean, double sd)
{
	std::normal_distribution<> d(mean,sd);
	return d(rng);
}

void RandomNumbers::poisson(std::vector<int>& random_list, double mean)
{
	std::poisson_distribution<int> d(mean);
	for (size_t i(0); i<random_list.size(); ++i)
	{
		random_list[i] = d(rng);
	}
}

int RandomNumbers::poisson(double mean)
{
	std::poisson_distribution<int> d(mean);
	return d(rng);
}



