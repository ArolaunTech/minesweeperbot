#include <vector>
#include <string>

#ifndef BOARD_H
#define BOARD_H

struct Board {
	std::vector<std::vector<int> > mines;

	bool isMine(int row, int col);
	int numMinesAround(int row, int col);

	void setMines(int rows, int cols, int nummines);

	std::string toString();
};

#endif