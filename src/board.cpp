#include <string>

#include "board.h"
#include "random.h"

void Board::calcMinesAround() {
	std::size_t numrows = mines.size();
	std::size_t numcols = mines[0].size();

	minesaround = std::vector<std::vector<int> >(numrows, std::vector<int>(numcols));

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (!isMine(i, j)) continue;

			for (int dr = -1; dr <= 1; dr++) {
				for (int dc = -1; dc <= 1; dc++) {
					if (dr == 0 && dc == 0) continue;

					int newi = i + dr;
					int newj = j + dc;

					if (newi < 0 || newi >= numrows || newj < 0 || newj >= numcols) continue;

					minesaround[newi][newj]++;
				}
			}
		}
	}
}

bool Board::isMine(int row, int col) {
	return mines[row][col] == 1;
}

int Board::numMinesAround(int row, int col) {
	return minesaround[row][col];
}

void Board::setMines(int rows, int cols, int nummines) {
	mines = std::vector<std::vector<int> >(rows, std::vector<int>(cols));

	for (int i = 0; i < nummines; i++) {
		int randrow, randcol;

		do {
			randrow = randint(0, rows - 1);
			randcol = randint(0, cols - 1);
		} while (mines[randrow][randcol]);

		mines[randrow][randcol] = 1;
	}

	calcMinesAround();
}

void Board::setMines(std::vector<std::vector<int> > newmines) {
	mines = newmines;

	calcMinesAround();
}

void Board::setMineStatus(int row, int col, bool minestatus) {
	mines[row][col] = minestatus ? 1 : 0;

	calcMinesAround();
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