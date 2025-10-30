#include <vector>
#include <set>

#include "game.h"

#ifndef SOLVER_H
#define SOLVER_H

struct TranspositionHash {
	std::size_t operator()(const std::pair<std::vector<int>, std::vector<bool> >& p) const {
		std::size_t h1 = 0;
        for (int x : p.first) {
            h1 ^= std::hash<int>{}(x) + 0x9e3779b9 + (h1 << 6) + (h1 >> 7);
        }

        std::size_t h2 = 0;
        for (bool b : p.second) {
            h2 ^= std::hash<bool>{}(b) + 0x7834dccb + (h2 << 10) + (h2 >> 1);
        }

        return h1 ^ (h2 << 2);
	}
};

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

struct SearchResult {
	std::size_t bestClick;
	int wins;
};

struct GameTree {
	BoardPosition move;

	std::vector<GameTree> children;
};

class Solver {
	std::vector<Move> queue;

	SearchResult search(
		Game& game, 
		std::vector<int> allowedIndices,
		std::vector<bool> allowedClicks,
		int totalCombinations, 
		std::vector<BoardPosition>& hiddenCells, 
		std::vector<std::vector<std::size_t> >& hiddenCellIndex,
		std::vector<std::vector<bool> >& mineCombinations,
		std::unordered_map<std::pair<std::vector<int>, std::vector<bool> >, SearchResult, TranspositionHash>& transpositionTable
	);
public:
	int bruteForceCalls;

	AnalyzeResult analyze(Game& game);
	PossibilitiesResult getPossibilities(Game& game);

	BoardPosition runBruteForce(Game& game);
	BoardPosition runGS(Game& game);
	bool runHeavyRollout(Game game);

	std::vector<std::vector<double> > getMineProbabilities(Game& game);
	Move getBestMove(Game& game);

	void clear();
};

#endif