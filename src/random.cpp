#include <random>

#include "random.h"

std::mt19937 generator(time(NULL));
std::uniform_int_distribution<int> randomint(0, 1000000000);
std::uniform_real_distribution<> uniform(0, 1);

int randint(int a, int b) {
	return a + (randomint(generator)) % (b - a + 1);
}

double randfloat(double a, double b) {
	return uniform(generator) * (b - a) + a;
}

bool randchance(double x) {
	return uniform(generator) < x;
}