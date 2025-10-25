#include <vector>
#include <iostream>

#include "debug.h"

template <typename T>
void logVector(std::vector<T> a) {
	for (std::size_t i = 0; i < a.size(); i++) {
		std::cout << a[i] << " ";
	}
}

template<typename T>
void logVector(std::vector<std::vector<T> > a) {
	for (std::size_t i = 0; i < a.size(); i++) {
		logVector(a[i]);
		std::cout << "\n";
	}
}