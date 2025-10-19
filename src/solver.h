#include <vector>

#include "game.h"

#ifndef SOLVER_H
#define SOLVER_H

struct BoardPosition {
	int row;
	int col;
};

bool operator<(const BoardPosition& lhs, const BoardPosition& rhs);
bool operator==(const BoardPosition& lhs, const BoardPosition& rhs);
bool operator>(const BoardPosition& lhs, const BoardPosition& rhs);

struct Move {
	int row;
	int col;
	bool flag;
};

struct AnalyzeResult {
	std::vector<std::vector<double> > probabilities;
	std::vector<std::vector<double> > information;
};

class Solver {
	std::vector<Move> queue;
public:
	AnalyzeResult analyze(Game game);
	std::vector<std::vector<double> > getMineProbabilities(Game game);
	Move getBestMove(Game game);

	void clear();
};

#endif