#include <vector>
#include <iostream>
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <array>

#include "game.h"
#include "solver.h"
#include "random.h"

template <typename T>
void logVector(std::vector<T> a) {
	for (std::size_t i = 0; i < a.size(); i++) {
		std::cout << a[i] << " ";
	}
}

double logFactorial(int n) {
	double out = 0;
	for (int i = 2; i <= n; i++) {
		out += std::log((double)i);
	}
	return out;
}

double lognCr(int n, int r) {
	if (r < 0) return -1e100;
	if (r > n) return -1e100;

	return logFactorial(n) - logFactorial(r) - logFactorial(n - r);
}

double logAdd(double logA, double logB) {
	double diff = std::abs(logA - logB);
	double maxLog = std::max(logA, logB);

	return maxLog + std::log(1 + std::exp(-diff));
}

bool operator<(const BoardPosition& lhs, const BoardPosition& rhs) {
	if (lhs.row < rhs.row) return true;
	if (lhs.row > rhs.row) return false;
	return lhs.col < rhs.col;
}

bool operator==(const BoardPosition& lhs, const BoardPosition& rhs) {
	return (lhs.row == rhs.row) && (lhs.col == rhs.col);
}

bool operator>(const BoardPosition& lhs, const BoardPosition& rhs) {
	if (lhs.row > rhs.row) return true;
	if (lhs.row < rhs.row) return false;
	return lhs.col > rhs.col;
}

int getCellsOfTypeAroundAnother(CellState type, Game& game, int i, int j, int radius = 1) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	int out = 0;

	for (int dr = -radius; dr <= radius; dr++) {
		for (int dc = -radius; dc <= radius; dc++) {
			if (dr == 0 && dc == 0) continue;

			int newrow = i + dr;
			int newcol = j + dc;

			if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) continue;

			if (game.getCell(newrow, newcol) == type) out++;
		}
	}

	return out;
}

bool isUnrelevant(Game& game, int row, int col) {
	if (game.getCell(row, col) != CELL_HIDDEN) return false;

	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	bool neighborless = true;
	for (int dr = -1; dr <= 1; dr++) {
		for (int dc = -1; dc <= 1; dc++) {
			if (dr == 0 && dc == 0) continue;

			int newr = row + dr;
			int newc = col + dc;

			if (newr < 0 || newr >= numrows) continue;
			if (newc < 0 || newc >= numcols) continue;

			CellState other = game.getCell(newr, newc);

			neighborless &= (other == CELL_HIDDEN || other == CELL_FLAG);
		}
	}

	return neighborless;
}

int getNumUnrelevantCells(Game& game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	int numUnrelevantCells = 0;

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (isUnrelevant(game, i, j)) numUnrelevantCells++;
		}
	}

	return numUnrelevantCells;
}

PossibilitiesResult Solver::getPossibilities(Game& game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();
	int totalmines = game.getTotalMines();

	/*=== Count flags ===*/
	int numflags = 0;
	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_FLAG) numflags++;
		}
	}

	/*=== Create cell groups ===*/
	std::vector<std::set<BoardPosition> > relevantCellGroups;
	std::vector<std::vector<BoardPosition> > relevantCells;
	std::vector<std::vector<std::ptrdiff_t> > groupIndices(
		numrows,
		std::vector<std::ptrdiff_t>(
			numcols,
			-1
	));
	std::vector<std::vector<std::set<std::ptrdiff_t> > > revealedCellGroups(
		numrows, 
		std::vector<std::set<std::ptrdiff_t> >(
			numcols
	));

	int numUnrelevantCells = getNumUnrelevantCells(game);

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) != CELL_HIDDEN) continue;

			std::set<BoardPosition> neighbors;
			for (int dr = -1; dr <= 1; dr++) {
				for (int dc = -1; dc <= 1; dc++) {
					if (dr == 0 && dc == 0) continue;

					int newr = i + dr;
					int newc = j + dc;

					if (newr < 0 || newr >= numrows) continue;
					if (newc < 0 || newc >= numcols) continue;

					CellState other = game.getCell(newr, newc);

					if (other == CELL_HIDDEN || other == CELL_FLAG) continue;

					neighbors.insert(BoardPosition {(int)newr, (int)newc});
				}
			}

			if (neighbors.size() == 0) continue;

			if (std::find(relevantCellGroups.begin(), relevantCellGroups.end(), neighbors) == relevantCellGroups.end()) {
				relevantCellGroups.push_back(neighbors);
				relevantCells.push_back(std::vector<BoardPosition>{BoardPosition {(int)i, (int)j}});
			} else {
				std::ptrdiff_t idx = std::find(relevantCellGroups.begin(), relevantCellGroups.end(), neighbors) - relevantCellGroups.begin();

				relevantCells[idx].push_back(BoardPosition {(int)i, (int)j});
			}

			std::ptrdiff_t idx = std::find(relevantCellGroups.begin(), relevantCellGroups.end(), neighbors) - relevantCellGroups.begin();

			groupIndices[i][j] = idx;

			for (const auto& neighbor : neighbors) {
				revealedCellGroups[neighbor.row][neighbor.col].insert(idx);
			}
		}
	}

	std::size_t numRelevantCellGroups = relevantCellGroups.size();

	/*=== Get possibilities ===*/
	std::vector<std::vector<std::size_t> > possibilities = {{}};

	for (std::size_t i = 0; i < numRelevantCellGroups; i++) {

		std::set<BoardPosition> constraints = relevantCellGroups[i];
		std::vector<BoardPosition> cellGroup = relevantCells[i];

		std::size_t numCells = cellGroup.size();

		std::vector<std::vector<std::size_t> > newPossibilities;
		for (const auto& possibility : possibilities) {

			int minMinesNeeded = 0;
			int maxMinesNeeded = numCells;

			for (const BoardPosition& constraint : constraints) {
				int minmines = 0;
				int maxmines = 0;
				int minesneeded = 
					(int)game.getCell(constraint.row, constraint.col) - 
					getCellsOfTypeAroundAnother(CELL_FLAG, game, constraint.row, constraint.col);

				std::set<std::ptrdiff_t> groupIndices = revealedCellGroups[constraint.row][constraint.col];

				for (const std::ptrdiff_t& groupIndex : groupIndices) {
					if (groupIndex == i) continue;

					if (groupIndex < i) {
						minmines += possibility[groupIndex];
						maxmines += possibility[groupIndex];
					} else {
						maxmines += relevantCells[groupIndex].size();
					}
				}

				minMinesNeeded = std::max(minMinesNeeded, minesneeded - maxmines);
				maxMinesNeeded = std::min(maxMinesNeeded, minesneeded - minmines);
			}

			if (maxMinesNeeded < minMinesNeeded) continue;

			for (int newmines = minMinesNeeded; newmines <= maxMinesNeeded; newmines++) {
				std::vector<std::size_t> newpossibility = possibility;
				newpossibility.push_back(newmines);

				newPossibilities.push_back(newpossibility);
			}
		}

		possibilities = newPossibilities;
	}

	/*=== Probability ===*/
	PossibilitiesResult out;

	out.logCombinations.clear();
	out.valid.clear();

	std::size_t numPossibilities = possibilities.size();
	double logTotalCombinations = -1e100;

	for (std::size_t i = 0; i < numPossibilities; i++) {

		double logCombinations = 0;
		std::size_t irrelevantMines = totalmines - numflags;

		for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
			logCombinations += lognCr(relevantCells[j].size(), possibilities[i][j]);

			irrelevantMines -= possibilities[i][j];
		}

		out.logCombinations.push_back(-1e100);
		out.valid.push_back(false);

		if (irrelevantMines > numUnrelevantCells) continue;
		if (irrelevantMines < 0) continue;

		logCombinations += lognCr(numUnrelevantCells, irrelevantMines);
		logTotalCombinations = logAdd(logTotalCombinations, logCombinations);

		out.logCombinations[i] = logCombinations;
		out.valid[i] = true;
	}

	out.possibilities = possibilities;
	out.groups = relevantCellGroups;
	out.logTotalCombinations = logTotalCombinations;
	out.relevantCells = relevantCells;
	
	return out;
}

AnalyzeResult Solver::analyze(Game& game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();
	int totalmines = game.getTotalMines();

	/*=== Count flags ===*/
	int numflags = 0;
	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_FLAG) numflags++;
		}
	}

	int numUnrelevantCells = getNumUnrelevantCells(game);

	/*=== Create cell groups ===*/
	PossibilitiesResult possibilities = getPossibilities(game);

	std::size_t numRelevantCellGroups = possibilities.groups.size();

	/*=== Probability ===*/
	std::size_t numPossibilities = possibilities.possibilities.size();

	std::vector<double> avgMines(numRelevantCellGroups, 0);
	std::vector<bool> allOnes(numRelevantCellGroups, true);
	double irrelevantMineDensity = 0;
	for (std::size_t i = 0; i < numPossibilities; i++) {

		if (!possibilities.valid[i]) continue;

		std::size_t irrelevantMines = totalmines - numflags;

		for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
			irrelevantMines -= possibilities.possibilities[i][j];
		}
		
		double probability = std::exp(possibilities.logCombinations[i] - possibilities.logTotalCombinations);
		for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
			avgMines[j] += probability * possibilities.possibilities[i][j] / possibilities.relevantCells[j].size();

			if (possibilities.possibilities[i][j] < possibilities.relevantCells[j].size()) allOnes[j] = false;

			if (avgMines[j] > 1) avgMines[j] = 1;
		}

		irrelevantMineDensity += probability * irrelevantMines / numUnrelevantCells;

		if (irrelevantMineDensity > 1) irrelevantMineDensity = 1;
	}

	for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
		if (allOnes[j]) avgMines[j] = 1;
	}

	std::vector<std::vector<double> > out(numrows, std::vector<double>(numcols));

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_HIDDEN) out[i][j] = irrelevantMineDensity;
			if (game.getCell(i, j) == CELL_FLAG || game.getCell(i, j) == CELL_MINE) out[i][j] = 1;
		}
	}

	for (std::size_t i = 0; i < numRelevantCellGroups; i++) {
		for (std::size_t j = 0; j < possibilities.relevantCells[i].size(); j++) {
			out[possibilities.relevantCells[i][j].row][possibilities.relevantCells[i][j].col] = avgMines[i];
		}
	}

	AnalyzeResult analytics;

	analytics.probabilities = out;

	return analytics;
}

BoardPosition Solver::runMCTS(Game& game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();
	std::size_t totalmines = game.getTotalMines();

	PossibilitiesResult possibilities = getPossibilities(game);

	/*=== Probabilities of each possibility ===*/
	std::size_t numRelevantCellGroups = possibilities.groups.size();

	int numUnrelevantCells = getNumUnrelevantCells(game);

	/*=== Probability ===*/
	std::size_t numPossibilities = possibilities.possibilities.size();
	std::vector<double> possibleProbs(numPossibilities, 0);

	for (std::size_t i = 0; i < numPossibilities; i++) {
		if (!possibilities.valid[i]) continue;
		
		possibleProbs[i] = std::exp(possibilities.logCombinations[i] - possibilities.logTotalCombinations);
	}

	/*=== Search ===*/
	for (int i = 0; i < 1000; i++) {
		double selector = randfloat(0, 1);
		int index = 0;

		while (selector > 0) {
			selector -= possibleProbs[index];
			index++;
		}

		index--;

		/*=== Generate mine possibility ===*/
		std::vector<std::size_t> possibility = possibilities.possibilities[index];
		std::vector<std::vector<bool> > mines(numrows, std::vector<bool>(numcols));

		for (std::size_t row = 0; row < numrows; row++) {
			for (std::size_t col = 0; col < numcols; col++) {
				mines[row][col] = game.getCell(row, col) == CELL_FLAG;
			}
		}

		for (std::size_t groupIndex = 0; groupIndex < numRelevantCellGroups; groupIndex++) {
			std::size_t groupmines = possibility[groupIndex];
			std::size_t assigned = 0;
			std::size_t numcells = possibilities.relevantCells[groupIndex].size();

			for (std::size_t cellIndex = 0; cellIndex < numcells; cellIndex++) {
				if (assigned == groupmines) break;

				bool mine = randfloat(0, 1) <= (double)(groupmines - assigned) / (numcells - cellIndex);

				if (mine) assigned++;

				mines[possibilities.relevantCells[groupIndex][cellIndex].row][possibilities.relevantCells[groupIndex][cellIndex].col] = mine;
			}
		}

		std::size_t irrelevantmines = totalmines;
		for (std::size_t row = 0; row < numrows; row++) {
			for (std::size_t col = 0; col < numcols; col++) {
				if (mines[row][col]) irrelevantmines--;
			}
		}

		std::size_t assigned = 0;
		std::size_t encountered = 0;
		for (std::size_t row = 0; row < numrows; row++) {
			for (std::size_t col = 0; col < numcols; col++) {
				if (!isUnrelevant(game, row, col)) continue;

				bool mine = randfloat(0, 1) <= (double)(irrelevantmines - assigned) / (numUnrelevantCells - encountered);

				if (mine) assigned++;
				encountered++;

				mines[row][col] = mine;
			}
		}

		/*=== Setup game ===*/
		std::vector<std::vector<bool> > revealed = game.getRevealed();
		std::vector<std::vector<bool> > flagged = game.getFlagged();

		logVector(possibility);
		std::cout << " mcts poss\n";

		for (std::size_t row = 0; row < numrows; row++) {
			for (std::size_t col = 0; col < numcols; col++) {
				std::cout << mines[row][col] << " ";
			}
			std::cout << "\n";
		}
	}

	return BoardPosition {0, 0};
}

std::vector<std::vector<double> > Solver::getMineProbabilities(Game& game) {
	return analyze(game).probabilities;
}

Move Solver::getBestMove(Game& game) {
	if (queue.size() > 0) {
		Move out = queue[queue.size() - 1];
		
		queue.pop_back();
		return out;
	}

	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	bool allHidden = true;

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			CellState cell = game.getCell(i, j);
			if (cell != CELL_HIDDEN) allHidden = false;

			if (cell == CELL_HIDDEN || cell == CELL_FLAG) continue;

			unsigned int numUnrevealedCells = getCellsOfTypeAroundAnother(CELL_HIDDEN, game, i, j);
			unsigned int numFlags = getCellsOfTypeAroundAnother(CELL_FLAG, game, i, j);

			if (numFlags == cell || numFlags + numUnrevealedCells == cell) {
				for (int dr = -1; dr <= 1; dr++) {
					for (int dc = -1; dc <= 1; dc++) {
						if (dr == 0 && dc == 0) continue;

						int newrow = i + dr;
						int newcol = j + dc;

						if (newrow < 0 || newrow >= numrows) continue;
						if (newcol < 0 || newcol >= numcols) continue;

						CellState other = game.getCell(newrow, newcol);

						if (other != CELL_HIDDEN) continue;

						queue.push_back(Move{newrow, newcol, (numFlags + numUnrevealedCells == cell)});
					}
				}
			}
		}
	}

	if (queue.size() > 0) {
		Move out = queue[queue.size() - 1];
		
		queue.pop_back();
		return out;
	}

	if (allHidden) {
		return Move{0, 0, false};
	}

	AnalyzeResult analytics = analyze(game);
	std::vector<std::vector<double> > probabilities = analytics.probabilities;

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) != CELL_HIDDEN) continue;

			if (probabilities[i][j] == 0) {
				queue.push_back(Move{(int)i, (int)j, false});
			} else if (probabilities[i][j] == 1) {
				queue.push_back(Move{(int)i, (int)j, true});
			}
		}
	}

	if (queue.size() > 0) {
		Move out = queue[queue.size() - 1];
		
		queue.pop_back();
		return out;
	}

	//BoardPosition bestguess = runMCTS(game);

	Move out;
	out.flag = false;

	double bestscore = -1;
	std::size_t bestrow = 0, bestcol = 0;
	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) != CELL_HIDDEN) continue;

			//Naive Rule 1 (37%)
			double score = -probabilities[i][j];

			if (score >= bestscore) {
				bestscore = score;
				bestrow = i;
				bestcol = j;
			}
		}
	}

	out.row = bestrow;
	out.col = bestcol;

	return out;
}

void Solver::clear() {
	queue.clear();
}