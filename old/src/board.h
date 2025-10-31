#include <vector>
#include <string>

#ifndef BOARD_H
#define BOARD_H

struct Board {
private:
	std::vector<std::vector<bool> > mines;
public:
	bool isMine(int row, int col);
	int numMinesAround(int row, int col);

	void setMines(int rows, int cols, int nummines);
	void setMines(std::vector<std::vector<bool> > newmines);
	void setMineStatus(int row, int col, bool minestatus);

	std::string toString();
};

#endif