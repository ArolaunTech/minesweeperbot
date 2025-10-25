#include <vector>
#include <string>

#include "board.h"

#ifndef GAME_H
#define GAME_H

enum CellState {
	CELL_HIDDEN = 9,
	CELL_MINE = 10,
	CELL_FLAG = 11,
	CELL_0MINES = 0,
	CELL_1MINES = 1,
	CELL_2MINES = 2,
	CELL_3MINES = 3,
	CELL_4MINES = 4,
	CELL_5MINES = 5,
	CELL_6MINES = 6,
	CELL_7MINES = 7,
	CELL_8MINES = 8,
};

enum GameState {
	GAME_WIN,
	GAME_LOSS,
	GAME_ONGOING
};

class Game {
	Board board;
	std::vector<std::vector<bool> > revealed;
	std::vector<std::vector<bool> > flagged;
	GameState state;
	int totalmines;

public:
	Game(int rows, int cols, int nummines);

	void initializeWithMineArrangement(std::vector<std::vector<bool> > mines);

	std::size_t getRows();
	std::size_t getCols();

	CellState getCell(int row, int col);
	bool cellIsState(int row, int col, CellState state);
	std::vector<std::vector<CellState> > getBoard();
	std::vector<std::vector<bool> > getRevealed();
	std::vector<std::vector<bool> > getFlagged();
	int getTotalMines();

	GameState getGameState();
	std::string toString();

	GameState click(int row, int col);
	void flag(int row, int col);
	void unflag(int row, int col);
};

#endif