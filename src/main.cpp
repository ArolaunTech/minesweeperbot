#include <iostream>
#include <string>

#include "board.h"

int main() {
	Board board;
	board.setMines(16, 30, 99);

	std::cout << board.toString();
}