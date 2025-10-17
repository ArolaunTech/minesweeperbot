#include <vector>
#include <iostream>
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>

#include "game.h"
#include "solver.h"

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

std::vector<std::vector<double> > Solver::getMineProbabilities(Game game) {
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
	std::vector<std::vector<std::set<std::ptrdiff_t> > > revealedCellGroups(
		numrows, 
		std::vector<std::set<std::ptrdiff_t> >(
			numcols
	));

	int numUnrelevantCells = 0;

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			CellState cell = game.getCell(i, j);

			if (cell != CELL_HIDDEN) continue;

			std::set<BoardPosition> neighbors;
			for (int dr = -1; dr <= 1; dr++) {
				for (int dc = -1; dc <= 1; dc++) {
					if (dr == 0 && dc == 0) continue;

					int newr = i + dr;
					int newc = j + dc;

					if (newr < 0 || newr >= numrows) continue;
					if (newc < 0 || newc >= numcols) continue;

					CellState other = game.getCell(newr, newc);

					if (other == CELL_HIDDEN || other == CELL_FLAG || other == CELL_MINE) continue;

					neighbors.insert(BoardPosition {(int)newr, (int)newc});
				}
			}

			if (neighbors.size() == 0) {
				numUnrelevantCells++;
				continue;
			}

			if (std::find(relevantCellGroups.begin(), relevantCellGroups.end(), neighbors) == relevantCellGroups.end()) {
				relevantCellGroups.push_back(neighbors);
				relevantCells.push_back(std::vector<BoardPosition>{BoardPosition {(int)i, (int)j}});
			} else {
				std::ptrdiff_t idx = std::find(relevantCellGroups.begin(), relevantCellGroups.end(), neighbors) - relevantCellGroups.begin();

				relevantCells[idx].push_back(BoardPosition {(int)i, (int)j});
			}

			for (const auto& neighbor : neighbors) {
				std::ptrdiff_t idx = std::find(relevantCellGroups.begin(), relevantCellGroups.end(), neighbors) - relevantCellGroups.begin();

				revealedCellGroups[neighbor.row][neighbor.col].insert(idx);
			}
		}
	}

	std::size_t numRelevantCellGroups = relevantCellGroups.size();

	/*=== No possibilities ===*/
	if (numRelevantCellGroups == 0) {
		std::vector<std::vector<double> > out(numrows, std::vector<double>(numcols, 0));

		for (std::size_t i = 0; i < numrows; i++) {
			for (std::size_t j = 0; j < numcols; j++) {
				if (game.getCell(i, j) == CELL_HIDDEN) {
					out[i][j] = (double)(totalmines - numflags) / numUnrelevantCells;
				} else if (game.getCell(i, j) == CELL_FLAG) {
					out[i][j] = 1;
				}
			}
		}

		return out;
	}

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
				int minesneeded = (int)game.getCell(constraint.row, constraint.col);

				for (int dr = -1; dr <= 1; dr++) {
					for (int dc = -1; dc <= 1; dc++) {
						if (dr == 0 && dc == 0) continue;	

						int newr = constraint.row + dr;
						int newc = constraint.col + dc;

						if (newr < 0 || newr >= numrows) continue;
						if (newc < 0 || newc >= numcols) continue;

						CellState other = game.getCell(newr, newc);

						if (other == CELL_FLAG) minesneeded--;
					}
				}

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
	std::size_t numPossibilities = possibilities.size();
	double logTotalCombinations = -1e100;

	for (std::size_t i = 0; i < numPossibilities; i++) {

		double logCombinations = 0;
		std::size_t irrelevantMines = totalmines - numflags;

		for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
			logCombinations += lognCr(relevantCells[j].size(), possibilities[i][j]);

			irrelevantMines -= possibilities[i][j];
		}

		if (irrelevantMines > numUnrelevantCells) continue;
		if (irrelevantMines < 0) continue;

		logCombinations += lognCr(numUnrelevantCells, irrelevantMines);
		logTotalCombinations = logAdd(logTotalCombinations, logCombinations);
	}

	std::vector<double> avgMines(numRelevantCellGroups, 0);
	double irrelevantMineDensity = 0;
	for (std::size_t i = 0; i < numPossibilities; i++) {

		double logCombinations = 0;
		std::size_t irrelevantMines = totalmines - numflags;

		for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
			logCombinations += lognCr(relevantCells[j].size(), possibilities[i][j]);

			irrelevantMines -= possibilities[i][j];
		}

		if (irrelevantMines > numUnrelevantCells) continue;
		if (irrelevantMines < 0) continue;

		logCombinations += lognCr(numUnrelevantCells, irrelevantMines);
		
		double probability = std::exp(logCombinations - logTotalCombinations);
		for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
			avgMines[j] += probability * possibilities[i][j] / relevantCells[j].size();

			if (avgMines[j] > 1) avgMines[j] = 1;
		}

		irrelevantMineDensity += probability * irrelevantMines / numUnrelevantCells;

		if (irrelevantMineDensity > 1) irrelevantMineDensity = 1;
	}

	std::vector<std::vector<double> > out(numrows, std::vector<double>(numcols));

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_HIDDEN) out[i][j] = irrelevantMineDensity;
			if (game.getCell(i, j) == CELL_FLAG || game.getCell(i, j) == CELL_MINE) out[i][j] = 1;
		}
	}

	for (std::size_t i = 0; i < numRelevantCellGroups; i++) {
		for (std::size_t j = 0; j < relevantCells[i].size(); j++) {
			out[relevantCells[i][j].row][relevantCells[i][j].col] = avgMines[i];
		}
	}

	return out;
}

/*
std::vector<std::vector<double> > Solver::getMineProbabilities(Game game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();
	int totalmines = game.getTotalMines();

	std::vector<std::vector<double> > out(numrows, std::vector<double>(numcols, 0));

	std::vector<std::vector<std::size_t> > unrevealedCellIdxs(numrows, std::vector<std::size_t>(numcols, numrows * numcols));
	std::vector<BoardPosition> relevantUnrevealedCells;
	int unrelevantCells = 0;

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_HIDDEN || game.getCell(i, j) == CELL_FLAG) {
				unrelevantCells++;
				continue;
			}
			if (game.getCell(i, j) == CELL_MINE) continue;

			for (int di = -1; di <= 1; di++) {
				for (int dj = -1; dj <= 1; dj++) {
					if (di == 0 && dj == 0) continue;

					int newi = i + di;
					int newj = j + dj;

					if (newi < 0 || newi >= numrows) continue;
					if (newj < 0 || newj >= numcols) continue;

					CellState other = game.getCell(newi, newj);

					if (other == CELL_HIDDEN || other == CELL_FLAG) {
						if (unrevealedCellIdxs[newi][newj] != numrows * numcols) continue;

						unrelevantCells--;

						unrevealedCellIdxs[newi][newj] = relevantUnrevealedCells.size();
						relevantUnrevealedCells.push_back(BoardPosition{(int)newi, (int)newj});
					}
				}
			}
		}
	}

	//std::cout << relevantUnrevealedCells.size() << "\n";

	if (relevantUnrevealedCells.empty()) {
		//The game has not started yet, each cell has an equal probability of a mine
		return std::vector<std::vector<double> >(numrows, std::vector<double>(numcols, (double)totalmines / (numrows * numcols)));
	}

	std::vector<std::vector<int> > possibilities = {{}};
	std::size_t numRelevantUnrevealedCells = relevantUnrevealedCells.size();

	for (std::size_t i = 0; i < numRelevantUnrevealedCells; i++) {
		int currrow = relevantUnrevealedCells[i].row;
		int currcol = relevantUnrevealedCells[i].col;

		std::size_t numpossibilities = possibilities.size();
		if (game.getCell(currrow, currcol) == CELL_FLAG) {
			for (std::size_t j = 0; j < numpossibilities; j++) {
				possibilities[j].push_back(1);
			}
			continue;
		}

		std::vector<std::vector<int> > newpossibilities;

		for (std::size_t j = 0; j < numpossibilities; j++) {
			for (int newmine = 0; newmine <= 1; newmine++) {
				std::vector<int> newpossibility = possibilities[j];
				newpossibility.push_back(newmine);

				//logVector(newpossibility);
				//std::cout << " poss\n\n";

				//std::cout << currrow << " " << currcol << "\n";

				bool valid = true;
				for (int dr = -1; dr <= 1; dr++) {
					for (int dc = -1; dc <= 1; dc++) {
						int newrow = currrow + dr, newcol = currcol + dc;

						if (newrow < 0 || newrow >= numrows) continue;
						if (newcol < 0 || newcol >= numcols) continue;

						CellState other = game.getCell(newrow, newcol);

						if (other == CELL_HIDDEN || other == CELL_MINE || other == CELL_FLAG) continue;

						int minmines = 0, maxmines = 0;

						for (int dr2 = -1; dr2 <= 1; dr2++) {
							for (int dc2 = -1; dc2 <= 1; dc2++) {
								int newrow2 = newrow + dr2;
								int newcol2 = newcol + dc2;

								if (newrow2 < 0 || newrow2 >= numrows) continue;
								if (newcol2 < 0 || newcol2 >= numcols) continue;

								CellState otherUnrevealed = game.getCell(newrow2, newcol2);

								if (otherUnrevealed == CELL_FLAG) {
									minmines++;
									maxmines++;
									continue;
								}
								if (otherUnrevealed != CELL_HIDDEN) continue;

								if (unrevealedCellIdxs[newrow2][newcol2] > i) {
									maxmines++;
								} else {
									minmines += newpossibility[unrevealedCellIdxs[newrow2][newcol2]];
									maxmines += newpossibility[unrevealedCellIdxs[newrow2][newcol2]];
								}
							}
						}

						//std::cout << newrow << " " << newcol << " " << minmines << " " << maxmines << "\n";

						if (other < minmines || other > maxmines) valid = false;
					}
				}

				if (valid) newpossibilities.push_back(newpossibility);
			}
		}

		possibilities = newpossibilities;

		//std::cout << possibilities.size() << "\n";
	}

	std::size_t numpossibilities = possibilities.size();
	double logTotalCombinations = -1e100;
	for (std::size_t i = 0; i < numpossibilities; i++) {
		int accountedMines = 0;
		for (std::size_t j = 0; j < numRelevantUnrevealedCells; j++) {
			accountedMines += possibilities[i][j];
		}

		if (accountedMines > totalmines) continue;
		if (totalmines - accountedMines > unrelevantCells) continue;

		logTotalCombinations = logAdd(logTotalCombinations, lognCr(unrelevantCells, totalmines - accountedMines));
	}

	//double logTotalGames = lognCr(numrows * numcols, totalmines);
	//std::cout << (logTotalGames - logTotalCombinations) / std::log(2) << " bits encoded in cell positions out of " << logTotalGames / std::log(2) << " total\n";

	std::vector<double> probabilities(numRelevantUnrevealedCells, 0);
	std::vector<bool> allOnes(numRelevantUnrevealedCells, true);
	double unrelevantProbability = 0;

	for (std::size_t i = 0; i < numpossibilities; i++) {
		//logVector(possibilities[i]);

		int accountedMines = 0;
		for (std::size_t j = 0; j < numRelevantUnrevealedCells; j++) {
			accountedMines += possibilities[i][j];
		}

		if (accountedMines > totalmines) continue;
		if (totalmines - accountedMines > unrelevantCells) continue;

		double probability = std::exp(lognCr(unrelevantCells, totalmines - accountedMines) - logTotalCombinations);

		for (std::size_t j = 0; j < numRelevantUnrevealedCells; j++) {
			if (possibilities[i][j] == 0) allOnes[j] = false;

			probabilities[j] += possibilities[i][j] * probability;
		}

		unrelevantProbability += probability * (totalmines - accountedMines) / unrelevantCells;

		//std::cout << " possfinal\n";
	}

	for (std::size_t j = 0; j < numRelevantUnrevealedCells; j++) {
		if (allOnes[j]) probabilities[j] = 1;
	}

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_FLAG) {
				out[i][j] = 1;
				continue;
			}
			if (game.getCell(i, j) != CELL_HIDDEN) continue;

			if (unrevealedCellIdxs[i][j] == numcols * numrows) {
				out[i][j] = unrelevantProbability;
				continue;
			}

			out[i][j] = probabilities[unrevealedCellIdxs[i][j]];
		}
	}

	return out;
}*/

Move Solver::getBestMove(Game game) {
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

			unsigned int numUnrevealedCells = 0;
			unsigned int numFlags = 0;

			for (int dr = -1; dr <= 1; dr++) {
				for (int dc = -1; dc <= 1; dc++) {
					if (dr == 0 && dc == 0) continue;

					int newrow = i + dr;
					int newcol = j + dc;

					if (newrow < 0 || newrow >= numrows) continue;
					if (newcol < 0 || newcol >= numcols) continue;

					CellState other = game.getCell(newrow, newcol);

					if (other == CELL_HIDDEN) numUnrevealedCells++;
					if (other == CELL_FLAG) numFlags++;
				}
			}

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

	std::vector<std::vector<double> > probabilities = getMineProbabilities(game);

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

	Move out;
	out.flag = false;

	//Naive rule #1 - always click on lowest probability (32% win rate)
	double lowestprobability = 1;
	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) != CELL_HIDDEN) continue;
			if (probabilities[i][j] < lowestprobability) lowestprobability = probabilities[i][j];
		}
	}

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (probabilities[i][j] > lowestprobability) continue;
			if (game.getCell(i, j) != CELL_HIDDEN) continue;

			out.row = i;
			out.col = j;
			return out;
		}
	}

	return out;
}

void Solver::clear() {
	queue.clear();
}