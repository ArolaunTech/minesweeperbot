#include <random>

#ifndef RANDOM_H
#define RANDOM_H

std::mt19937 generator(time(NULL));
std::uniform_int_distribution<int> randomint(0, 1000000000);

int randint(int a, int b) {
	return a + (randomint(generator)) % (b - a + 1);
}

#endif