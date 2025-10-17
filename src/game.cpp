#include <vector>

#include "game.h"
#include "board.h"
#include "random.h"

Game::Game(int rows, int cols, int nummines) {
	board.setMines(rows, cols, nummines);

	revealed = std::vector<std::vector<bool> >(rows, std::vector<bool>(cols, false));
	flagged = std::vector<std::vector<bool> >(rows, std::vector<bool>(cols, false));

	state = GAME_ONGOING;

	totalmines = nummines;
}

std::size_t Game::getRows() {
	return revealed.size();
}

std::size_t Game::getCols() {
	return revealed[0].size();
}

CellState Game::getCell(int row, int col) {
	if (flagged[row][col]) return CELL_FLAG;
	if (!revealed[row][col]) return CELL_HIDDEN;
	if (board.isMine(row, col)) return CELL_MINE;

	return static_cast<CellState>(board.numMinesAround(row, col));
}

std::vector<std::vector<CellState> > Game::getBoard() {
	std::size_t numrows = getRows();
	std::size_t numcols = getCols();

	std::vector<std::vector<CellState> > out(numrows, std::vector<CellState>(numcols));

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (!revealed[i][j]) out[i][j] = CELL_HIDDEN;
			else if (board.isMine(i, j)) out[i][j] = CELL_MINE;
			else out[i][j] = static_cast<CellState>(board.numMinesAround(i, j));
		}
	}

	return out;
}

std::vector<std::vector<bool> > Game::getRevealed() {
	return revealed;
}

std::vector<std::vector<bool> > Game::getFlagged() {
	return flagged;
}

int Game::getTotalMines() {
	return totalmines;
}

GameState Game::getGameState() {
	return state;
}

std::string Game::toString() {
	std::string out;

	switch (state) {
	case GAME_WIN:
		out = "Game won\n";
		break;
	case GAME_LOSS:
		out = "Game lost\n";
		break;
	case GAME_ONGOING:
		out = "Game ongoing\n";
		break;
	}

	std::size_t numrows = getRows();
	std::size_t numcols = getCols();

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (revealed[i][j]) {
				int minesaround = board.numMinesAround(i, j);

				if (board.isMine(i,j) == 1) {
					out += "M";
				} else if (minesaround > 0) {
					out += std::to_string(minesaround);
				} else {
					out += " ";
				}
			} else if (flagged[i][j]) {
				out += "F";
			} else {
				out += ".";
			}

			out += " ";
		}
		out += "\n";
	}

	return out;
}



GameState Game::click(int row, int col) {
	if (flagged[row][col]) return state;
	if (revealed[row][col]) return state;

	std::size_t numrows = getRows();
	std::size_t numcols = getCols();

	if (board.isMine(row, col)) {
		for (std::size_t i = 0; i < numrows; i++) {
			for (std::size_t j = 0; j < numcols; j++) {
				if (revealed[i][j]) {
					revealed[row][col] = true;
					state = GAME_LOSS;
					return state;
				}
			}
		}

		int randrow, randcol;

		do {
			randrow = randint(0, numrows - 1);
			randcol = randint(0, numcols - 1);
		} while (board.isMine(randrow, randcol));

		board.setMineStatus(randrow, randcol, true);
		board.setMineStatus(row, col, false);

		return click(row, col);
	}

	revealed[row][col] = true;

	int minesaround = board.numMinesAround(row, col);
	if (minesaround == 0) {
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				int newrow = row + i;
				int newcol = col + j;
				if (newrow < 0 || newrow >= numrows) continue;
				if (newcol < 0 || newcol >= numcols) continue;

				click(newrow, newcol);
			}
		}
	}

	for (int i = 0; i < numrows; i++) {
		for (int j = 0; j < numcols; j++) {
			if (!revealed[i][j] && !board.isMine(i, j)) return state;
		}
	}
	state = GAME_WIN;
	return state;
}

void Game::flag(int row, int col) {
	flagged[row][col] = true;
}

void Game::unflag(int row, int col) {
	flagged[row][col] = false;
}