#include <vector>
#include <iostream>
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <array>

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

int getCellsOfTypeAroundAnother(CellState type, Game game, int i, int j) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	int out = 0;

	for (int dr = -1; dr <= 1; dr++) {
		for (int dc = -1; dc <= 1; dc++) {
			if (dr == 0 && dc == 0) continue;

			int newrow = i + dr;
			int newcol = j + dc;

			if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) continue;

			if (game.getCell(newrow, newcol) == type) out++;
		}
	}

	return out;
}

AnalyzeResult Solver::analyze(Game game) {
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
	std::vector<bool> allOnes(numRelevantCellGroups, true);
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

			if (possibilities[i][j] < relevantCells[j].size()) allOnes[j] = false;

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
		for (std::size_t j = 0; j < relevantCells[i].size(); j++) {
			out[relevantCells[i][j].row][relevantCells[i][j].col] = avgMines[i];
		}
	}

	/*=== Get average info returns ===*/
	std::vector<std::vector<double> > avgInfoReturns(numrows, std::vector<double>(numcols, 0));

	//Get info returns for cells in middle of nowhere (optimization)
	double logIsolatedCombinations = -1e100;
	std::array<double, 9> logIsolatedResults;

	logIsolatedResults.fill(-1e100);

	for (std::size_t i = 0; i < numPossibilities; i++) {
		double logCombinations = 0;
		std::size_t irrelevantMines = totalmines - numflags;

		for (std::size_t j = 0; j < numRelevantCellGroups; j++) {
			logCombinations += lognCr(relevantCells[j].size(), possibilities[i][j]);

			irrelevantMines -= possibilities[i][j];
		}

		if (irrelevantMines >= numUnrelevantCells) continue;
		if (irrelevantMines < 0) continue;

		for (int nummines = 0; nummines <= 8; nummines++) {
			double logCombinationsFinal = logCombinations + lognCr(numUnrelevantCells - 9, irrelevantMines - nummines) + lognCr(8, nummines);

			logIsolatedResults[nummines] = logAdd(logIsolatedResults[nummines], logCombinationsFinal);
			logIsolatedCombinations = logAdd(logIsolatedCombinations, logCombinationsFinal);
		}
	}

	double isolatedInfoReturns = 0;
	for (int nummines = 0; nummines <= 8; nummines++) {
		isolatedInfoReturns += (logIsolatedCombinations - logIsolatedResults[nummines]) * std::exp(logIsolatedResults[nummines] - logIsolatedCombinations);
	}

	//Get all info returns
	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) != CELL_HIDDEN) continue;

			bool isolated = true;
			for (int dr = -2; dr <= 2; dr++) {
				for (int dc = -2; dc <= 2; dc++) {
					int newrow = i + dr;
					int newcol = j + dc;

					if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) continue;

					isolated &= game.getCell(newrow, newcol) == CELL_HIDDEN;
				}
			}
			if (isolated) {
				avgInfoReturns[i][j] = isolatedInfoReturns;
				continue;
			}

			std::ptrdiff_t currGroupIndex = groupIndices[i][j];

			int freeCells = 0;
			std::vector<int> neighboringGroups(numRelevantCellGroups + 1);
			for (int dr = -1; dr <= 1; dr++) {
				for (int dc = -1; dc <= 1; dc++) {
					if (dr == 0 && dc == 0) continue;

					int newrow = i + dr;
					int newcol = j + dc;

					if (newrow < 0 || newcol < 0 || newrow >= numrows || newcol >= numcols) continue;

					if (game.getCell(newrow, newcol) != CELL_HIDDEN) continue;

					freeCells++;

					std::ptrdiff_t groupIndex = groupIndices[newrow][newcol];

					if (groupIndex == -1) {
						neighboringGroups[numRelevantCellGroups]++;
					} else {
						neighboringGroups[groupIndex]++;
					}
				}
			}

			if (freeCells == 0) continue;

			double logTotalCombinations = -1e100;
			std::array<double, 9> logCombinations;
			logCombinations.fill(-1e100);

			for (std::size_t k = 0; k < numPossibilities; k++) {
				std::vector<int> minmines(numRelevantCellGroups + 1);
				std::vector<int> maxmines(numRelevantCellGroups + 1);

				bool valid = true;

				for (std::size_t groupIndex = 0; groupIndex < numRelevantCellGroups; groupIndex++) {
					std::size_t numCells = relevantCells[groupIndex].size();
					std::size_t numMines = possibilities[k][groupIndex];
					int numTouching = neighboringGroups[groupIndex];

					if (groupIndex == currGroupIndex) numCells--;
					if (numMines > numCells) {
						valid = false;
						break;
					}

					maxmines[groupIndex] = std::min(numTouching, (int)numMines);
					minmines[groupIndex] = std::max(0, (int)(numMines + numTouching - numCells));
				}

				if (!valid) continue;

				std::size_t irrelevantMines = totalmines - numflags;

				for (std::size_t groupIndex = 0; groupIndex < numRelevantCellGroups; groupIndex++) {
					irrelevantMines -= possibilities[k][groupIndex];
				}

				int numCells = numUnrelevantCells;
				if (currGroupIndex == -1) numCells--;

				if (irrelevantMines > numCells) continue;

				maxmines[numRelevantCellGroups] = std::min(neighboringGroups[numRelevantCellGroups], (int)irrelevantMines);
				minmines[numRelevantCellGroups] = std::max(0, (int)(irrelevantMines + neighboringGroups[numRelevantCellGroups] - numCells));

				/*=== Find combinations ===*/
				int numcombinations = 1;
				for (std::size_t groupIndex = 0; groupIndex <= numRelevantCellGroups; groupIndex++) {
					numcombinations *= maxmines[groupIndex] - minmines[groupIndex] + 1;
				}

				for (int combinationIndex = 0; combinationIndex < numcombinations; combinationIndex++) {
					double logPossibilityCombinations = 0;
					int cellNumber = 0;

					int curr = combinationIndex;
					for (std::size_t groupIndex = 0; groupIndex < numRelevantCellGroups; groupIndex++) {
						int groupMines = minmines[groupIndex] + curr % (maxmines[groupIndex] - minmines[groupIndex] + 1);
						std::size_t numCells = relevantCells[groupIndex].size();
						
						if (groupIndex == currGroupIndex) numCells--;					

						logPossibilityCombinations += 
							lognCr(neighboringGroups[groupIndex], groupMines) + 
							lognCr(numCells - neighboringGroups[groupIndex], possibilities[k][groupIndex] - groupMines);

						curr /= maxmines[groupIndex] - minmines[groupIndex] + 1;

						cellNumber += groupMines;
					}

					int groupMines = minmines[numRelevantCellGroups] + curr % (maxmines[numRelevantCellGroups] - minmines[numRelevantCellGroups] + 1);
					cellNumber += groupMines;

					int numCells = numUnrelevantCells;
					if (currGroupIndex == -1) numCells--;

					logPossibilityCombinations += 
						lognCr(neighboringGroups[numRelevantCellGroups], groupMines) + 
						lognCr(numCells - neighboringGroups[numRelevantCellGroups], irrelevantMines - groupMines);

					logCombinations[cellNumber] = logAdd(logCombinations[cellNumber], logPossibilityCombinations);

					//std::cout << logCombinations[cellNumber] << " " << logPossibilityCombinations << "\n";
				}
			}

			for (int cellNumber = 0; cellNumber <= 8; cellNumber++) {
				logTotalCombinations = logAdd(logTotalCombinations, logCombinations[cellNumber]);
			}

			double totalInfo = 0;

			for (int cellNumber = 0; cellNumber <= 8; cellNumber++) {
				totalInfo += (logTotalCombinations - logCombinations[cellNumber]) * std::exp(logCombinations[cellNumber] - logTotalCombinations);
			}

			avgInfoReturns[i][j] = totalInfo;
		}
	}

	AnalyzeResult analytics;

	analytics.probabilities = out;
	analytics.information = avgInfoReturns;

	return analytics;
}

std::vector<std::vector<double> > Solver::getMineProbabilities(Game game) {
	return analyze(game).probabilities;
}

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

	AnalyzeResult analytics = analyze(game);
	std::vector<std::vector<double> > probabilities = analytics.probabilities;
	std::vector<std::vector<double> > info = analytics.information;

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

	//Naive rule #2 - always click on highest (info / probability)
	double bestscore = -1;
	std::size_t bestrow = 0, bestcol = 0;
	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) != CELL_HIDDEN) continue;

			//Naive Rule 1 (37%)
			double score = -probabilities[i][j];

			//Naive Rule 2 (34%)
			//double score = info[i][j] / probabilities[i][j];

			if (score > bestscore) {
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