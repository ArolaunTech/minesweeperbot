#include <string>

#include "board.h"
#include "random.h"

bool Board::isMine(int row, int col) {
	return mines[row][col];
}

int Board::numMinesAround(int row, int col) {
	int out = 0;
	std::size_t numrows = mines.size();
	std::size_t numcols = mines[0].size();
	
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			if (i == 0 && j == 0) continue;

			int newrow = row + i;
			int newcol = col + j;
			if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) continue;

			out += mines[newrow][newcol] ? 1 : 0;
		}
	}

	return out;
}

void Board::setMines(int rows, int cols, int nummines) {
	mines = std::vector<std::vector<bool> >(rows, std::vector<bool>(cols, false));

	for (int i = 0; i < nummines; i++) {
		int randrow, randcol;

		do {
			randrow = randint(0, rows - 1);
			randcol = randint(0, cols - 1);
		} while (mines[randrow][randcol]);

		mines[randrow][randcol] = true;
	}
}

void Board::setMines(std::vector<std::vector<bool> > newmines) {
	mines = newmines;
}

void Board::setMineStatus(int row, int col, bool minestatus) {
	mines[row][col] = minestatus;
}

std::string Board::toString() {
	std::string out;

	std::size_t numrows = mines.size();
	std::size_t numcols = mines[0].size();

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			int minesaround = numMinesAround(i, j);

			if (mines[i][j] == 1) {
				out += "M";
			} else if (minesaround > 0) {
				out += std::to_string(minesaround);
			} else {
				out += " ";
			}

			out += " ";
		}
		out += "\n";
	}

	return out;
}