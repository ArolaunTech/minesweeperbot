#include "math.h"

int factorial(int n) {
	int out = 1;
	for (int i = 2; i <= n; i++) {
		out *= i;
	}
	return out;
}

int nCr(int n, int r) {
	if (r < 0 || r > n) return 0;

	int out = 1;

	for (int i = 1; i <= r; i++) {
		out *= n - r + i;
		out /= i;
	}

	return out;
}