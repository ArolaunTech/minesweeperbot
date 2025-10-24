#include <vector>
#include <set>

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

struct PossibilitiesResult {
	std::vector<std::set<BoardPosition> > groups;
	std::vector<std::vector<BoardPosition> > relevantCells;

	std::vector<std::vector<std::size_t> > possibilities;
	std::vector<bool> valid;
	std::vector<double> logCombinations;
	double logTotalCombinations;
};

struct AnalyzeResult {
	PossibilitiesResult possibilities;

	std::vector<std::vector<double> > probabilities;
};

class Solver {
	std::vector<Move> queue;
public:
	AnalyzeResult analyze(Game& game);
	PossibilitiesResult getPossibilities(Game& game);

	BoardPosition runBruteForce(Game& game);

	std::vector<std::vector<double> > getMineProbabilities(Game& game);
	Move getBestMove(Game& game);

	void clear();
};

#endif