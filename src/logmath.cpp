#include <cmath>
#include <algorithm>

#include "logmath.h"

double logFactorial(int n) {
	double out = 0;
	for (int i = 2; i <= n; i++) {
		out += std::log((double)i);
	}
	return out;
}

double lognCr(int n, int r) {
	if (r < 0) return -1e100;
	if (r > n) return -1e100;

	double out = 0;

	for (int i = n - r + 1; i <= n; i++) out += std::log((double)i);

	return out - logFactorial(r);
}

double logAdd(double logA, double logB) {
	double diff = std::abs(logA - logB);
	double maxLog = std::max(logA, logB);

	return maxLog + std::log(1 + std::exp(-diff));
}