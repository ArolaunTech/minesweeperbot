#include <string>

#include "board.h"
#include "random.h"

bool Board::isMine(int row, int col) {
	return mines[row][col] == 1;
}

int Board::numMinesAround(int row, int col) {
	int out = 0;
	std::size_t numrows = mines.size();
	std::size_t numcols = mines[0].size();
	
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			if (i == 0 && j == 0) {
				continue;
			}

			int newrow = row + i;
			int newcol = col + j;
			if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) {
				continue;
			}

			out += mines[newrow][newcol];
		}
	}

	return out;
}

void Board::setMines(int rows, int cols, int nummines) {
	mines.clear();

	for (int i = 0; i < rows; i++) {
		mines.push_back(std::vector<int>(cols));
	}

	for (int i = 0; i < nummines; i++) {
		int randrow, randcol;

		do {
			randrow = randint(0, rows - 1);
			randcol = randint(0, cols - 1);
		} while (mines[randrow][randcol] == 1);

		mines[randrow][randcol] = 1;
	}
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