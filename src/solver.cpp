#include <vector>
#include <iostream>
#include <cmath>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <array>
#include <utility>
#include <functional>

#include "game.h"
#include "solver.h"
#include "random.h"
#include "logmath.h"
#include "debug.h"

struct TranspositionHash {
	std::size_t operator()(const std::pair<std::vector<int>, std::vector<bool> >& p) const {
		std::size_t h1 = 0;
        for (int x : p.first) {
            h1 ^= std::hash<int>{}(x) + 0x9e3779b9 + (h1 << 6) + (h1 >> 7); // Hash combining for vector<int>
        }

        std::size_t h2 = 0;
        for (bool b : p.second) {
            h2 ^= std::hash<bool>{}(b) + 0x7834dccb + (h2 << 10) + (h2 >> 1); // Hash combining for vector<bool>
        }

        // Combine the two hashes
        return h1 ^ (h2 << 2); // Simple combination, more robust methods exist
	}
};

int factorial(int n) {
	int out = 1;
	for (int i = 2; i <= n; i++) {
		out *= i;
	}
	return out;
}

int nCr(int n, int r) {
	if (r < 0 || r > n) return 0;

	int out = 1;

	for (int i = 1; i <= r; i++) {
		out *= n - r + i;
		out /= i;
	}

	return out;
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

			if (game.cellIsState(newrow, newcol, type)) out++;
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

int getNumFlags(Game& game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	int numflags = 0;

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_FLAG) numflags++;
		}
	}

	return numflags;
}

int getNumHidden(Game& game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	int numhidden = 0;

	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (game.getCell(i, j) == CELL_HIDDEN) numhidden++;
		}
	}

	return numhidden;
}

bool Solver::runHeavyRollout(Game game) {
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	while (game.getGameState() == GAME_ONGOING) {
		bool moved = false;
		for (std::size_t i = 0; i < numrows; i++) {
			for (std::size_t j = 0; j < numcols; j++) {
				CellState cell = game.getCell(i, j);
				if (cell != CELL_HIDDEN) continue;

				int numUnrevealedCells = getCellsOfTypeAroundAnother(CELL_HIDDEN, game, i, j);
				int numFlags = getCellsOfTypeAroundAnother(CELL_FLAG, game, i, j);
				
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

							moved = true;

							if (numFlags + numUnrevealedCells == cell) {
								game.flag(newrow, newcol);
							} else {
								game.click(newrow, newcol);
							}
						}
					}
				}
			}
		}

		if (moved) continue;

		AnalyzeResult analytics = analyze(game);
		std::vector<std::vector<double> > probabilities = analytics.probabilities;

		double bestscore = 1;
		std::size_t bestrow = 0;
		std::size_t bestcol = 0;
		for (std::size_t i = 0; i < numrows; i++) {
			for (std::size_t j = 0; j < numcols; j++) {
				if (game.getCell(i, j) != CELL_HIDDEN) continue;

				if (probabilities[i][j] == 0) {
					game.click(i, j);
					moved = true;
				} else if (probabilities[i][j] == 1) {
					game.flag(i, j);
					moved = true;
				} else if (probabilities[i][j] < bestscore) {
					bestscore = probabilities[i][j];
					bestrow = i;
					bestcol = j;
				}
			}
		}

		if (moved) continue;

		game.click(bestrow, bestcol);
	}

	return game.getGameState() == GAME_WIN;
}

PossibilitiesResult Solver::getPossibilities(Game& game) {
	//TODO: Optimize. Currently this is the hottest function in the solver.
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();
	int totalmines = game.getTotalMines();

	/*=== Count flags ===*/
	int numflags = getNumFlags(game);

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
	std::vector<std::vector<std::size_t> > newPossibilities;

	for (std::size_t i = 0; i < numRelevantCellGroups; i++) {

		std::set<BoardPosition> constraints = relevantCellGroups[i];
		std::vector<BoardPosition> cellGroup = relevantCells[i];

		std::size_t numCells = cellGroup.size();

		newPossibilities.clear();
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

		possibilities.swap(newPossibilities);
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
	int numflags = getNumFlags(game);
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
	analytics.possibilities = possibilities;

	return analytics;
}

std::vector<std::vector<double> > Solver::getMineProbabilities(Game& game) {
	return analyze(game).probabilities;
}

std::vector<bool> getIthMineCombination(int numcells, int nummines, int i) {
	if (numcells <= 0) return std::vector<bool>();
	if (numcells == 1) return std::vector<bool>{nummines == 1};

	int combs0 = nCr(numcells - 1, nummines);

	bool currmine = i >= combs0;
	std::vector<bool> out = getIthMineCombination(numcells - 1, nummines - (currmine ? 1 : 0), i - (currmine ? combs0 : 0));

	out.push_back(currmine);

	return out;
}

struct BruteForceSearchResult {
	std::size_t bestClick;
	int wins;
};

BruteForceSearchResult bruteForceSearch(
	Game& game, 
	std::vector<int> allowedIndices,
	std::vector<bool> allowedClicks,
	int totalCombinations, 
	std::vector<BoardPosition>& hiddenCells, 
	std::vector<std::vector<std::size_t> >& hiddenCellIndex,
	std::vector<std::vector<bool> >& mineCombinations,
	std::unordered_map<std::pair<std::vector<int>, std::vector<bool> >, BruteForceSearchResult, TranspositionHash>& transpositionTable
) {
	if (transpositionTable.contains(std::pair<std::vector<int>, std::vector<bool> >(allowedIndices, allowedClicks))) {
		return transpositionTable[std::pair<std::vector<int>, std::vector<bool> >(allowedIndices, allowedClicks)];
	}

	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	std::size_t numHiddenCells = hiddenCells.size();

	std::size_t bestClick = 0;
	int bestWins = 0;
	bool forcedClick = false;

	for (std::size_t click = 0; click < numHiddenCells; click++) {
		if (!allowedClicks[click]) continue;

		int losses = 0;

		for (int i : allowedIndices) {
			if (mineCombinations[i][click]) losses++;
		}

		if (losses == 0) {
			bestClick = click;
			forcedClick = true;
			break;
		}
	}

	for (std::size_t click = 0; click < numHiddenCells; click++) {
		if (!allowedClicks[click]) continue;
		if (forcedClick && click != bestClick) continue;

		std::vector<bool> newAllowedClicks = allowedClicks;
		newAllowedClicks[click] = false;

		BoardPosition clickPosition = hiddenCells[click];

		std::vector<std::vector<int> > indices(9);
		int potentialwins = 0;

		for (int i : allowedIndices) {
			if (mineCombinations[i][click]) continue;

			int cellnumber = 0;
			potentialwins++;

			for (int dr = -1; dr <= 1; dr++) {
				for (int dc = -1; dc <= 1; dc++) {
					if (dr == 0 && dc == 0) continue;

					int newrow = clickPosition.row + dr;
					int newcol = clickPosition.col + dc;

					if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) continue;

					CellState other = game.getCell(newrow, newcol);

					if (other == CELL_FLAG) {
						cellnumber++;
						continue;
					}

					if (other == CELL_HIDDEN && mineCombinations[i][hiddenCellIndex[newrow][newcol]]) {
						cellnumber++;
					}
				}
			}

			indices[cellnumber].push_back(i);
		}

		if (potentialwins <= bestWins) continue;

		int wins = 0;
		for (int i = 0; i < 9; i++) {
			std::size_t numRemainingPossibilities = indices[i].size();

			if (numRemainingPossibilities == 0) continue;
			if (numRemainingPossibilities == 1) {
				wins++;
				continue;
			}

			wins += bruteForceSearch(
				game, 
				indices[i], 
				newAllowedClicks, 
				totalCombinations, 
				hiddenCells, 
				hiddenCellIndex, 
				mineCombinations,
				transpositionTable
			).wins;

			//std::cout << wins << " " << i << " wins\n\n";
		}
		//logVector(indices);
		//std::cout << "indices\n";
		//logVector(allowedClicks);
		//std::cout << " allowed\n";

		//std::cout << wins << " wins\n";
		//std::cout << forcedClick << " forced\n\n";

		if (wins > bestWins) {
			bestClick = click;
			bestWins = wins;
		}
	}

	BruteForceSearchResult out;
	out.bestClick = bestClick;
	out.wins = bestWins;

	transpositionTable[std::pair<std::vector<int>, std::vector<bool> >(allowedIndices, allowedClicks)] = out;

	return out;
}

BoardPosition Solver::runBruteForce(Game& game) {
	/*
	 * Consider all possible mine combinations and choose a perfect move based on that.
	 * Should only run near the end of the game when there are relatively few moves left.
	 */

	/*=== Get basic game info ===*/
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();
	int totalmines = game.getTotalMines();

	int numUnrelevantCells = getNumUnrelevantCells(game);
	int numflags = getNumFlags(game);

	std::vector<BoardPosition> hiddenCells;
	std::vector<std::vector<std::size_t> > hiddenCellIndex(numrows, std::vector<std::size_t>(numcols, -1));

	for (std::size_t row = 0; row < numrows; row++) {
		for (std::size_t col = 0; col < numcols; col++) {
			if (game.getCell(row, col) != CELL_HIDDEN) continue;

			hiddenCellIndex[row][col] = hiddenCells.size();
			hiddenCells.push_back(BoardPosition{(int)row, (int)col});
		}
	}

	std::size_t numHiddenCells = hiddenCells.size();

	/*=== Analyze game ===*/
	AnalyzeResult analytics = analyze(game);
	PossibilitiesResult possibilities = analytics.possibilities;

	/*=== Enumerate combinations ===*/
	std::size_t numPossibilities = possibilities.possibilities.size();
	std::size_t numRelevantCellGroups = possibilities.relevantCells.size();
	int totalCombinations = 0;

	std::vector<std::vector<bool> > mineCombinations;

	for (std::size_t i = 0; i < numPossibilities; i++) {
		if (!possibilities.valid[i]) continue;

		/*=== Calculate number of combinations for possibility ===*/
		int possibilityCombinations = 1;
		int irrelevantMines = totalmines - numflags;

		std::vector<int> groupDivs(numRelevantCellGroups + 1);
		std::vector<int> groupMods(numRelevantCellGroups);

		for (std::size_t groupIndex = 0; groupIndex < numRelevantCellGroups; groupIndex++) {
			int groupCombinations = nCr(
				possibilities.relevantCells[groupIndex].size(), 
				possibilities.possibilities[i][groupIndex]
			);

			groupMods[groupIndex] = groupCombinations;
			groupDivs[groupIndex] = possibilityCombinations;

			possibilityCombinations *= groupCombinations;

			irrelevantMines -= possibilities.possibilities[i][groupIndex];
		}

		groupDivs[numRelevantCellGroups] = possibilityCombinations;

		possibilityCombinations *= nCr(numUnrelevantCells, irrelevantMines);
		totalCombinations += possibilityCombinations;

		/*=== Store all possible combinations ===*/
		for (int possIndex = 0; possIndex < possibilityCombinations; possIndex++) {
			/*=== Initialize board ===*/
			std::vector<bool> mines(numHiddenCells);

			for (std::size_t groupIndex = 0; groupIndex < numRelevantCellGroups; groupIndex++) {
				int groupCombinationIndex = (possIndex / groupDivs[groupIndex]) % groupMods[groupIndex];

				std::vector<bool> mineCombo = getIthMineCombination(
					possibilities.relevantCells[groupIndex].size(), 
					possibilities.possibilities[i][groupIndex], 
					groupCombinationIndex
				);

				for (std::size_t cellIndex = 0; cellIndex < possibilities.relevantCells[groupIndex].size(); cellIndex++) {
					mines[
						hiddenCellIndex[
							possibilities.relevantCells[groupIndex][cellIndex].row
						][
							possibilities.relevantCells[groupIndex][cellIndex].col
						]
					] =
						mineCombo[cellIndex];
				}
			}

			int groupCombinationIndex = possIndex / groupDivs[numRelevantCellGroups];

			std::vector<bool> mineCombo = getIthMineCombination(
				numUnrelevantCells, 
				irrelevantMines, 
				groupCombinationIndex
			);

			int seen = 0;

			for (std::size_t row = 0; row < numrows; row++) {
				for (std::size_t col = 0; col < numcols; col++) {
					if (!isUnrelevant(game, row, col)) continue;

					mines[hiddenCellIndex[row][col]] = mineCombo[seen];

					seen++;
				}
			}

			mineCombinations.push_back(mines);
		}
	}

	//std::vector<std::vector<bool> >(1 << numHiddenCells, std::vector<bool>())

	std::vector<int> allowedIndices(totalCombinations);
	for (int i = 0; i < totalCombinations; i++) {
		allowedIndices.push_back(i);
	}

	std::unordered_map<std::pair<std::vector<int>, std::vector<bool> >, BruteForceSearchResult, TranspositionHash> transpositionTable;

	BruteForceSearchResult result = bruteForceSearch(
		game, 
		allowedIndices,
		std::vector<bool>(numHiddenCells, true),
		totalCombinations,
		hiddenCells,
		hiddenCellIndex,
		mineCombinations,
		transpositionTable
	);

	return hiddenCells[result.bestClick];
}

BoardPosition Solver::runGS(Game& game) {
	//Get best guessing move without brute force
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	AnalyzeResult analytics = analyze(game);

	std::vector<std::vector<double> > probabilities = analytics.probabilities;

	std::size_t bestrow = 0, bestcol = 0;



	/*=== Keep this for now ===*/
	double bestscore = -1;
	bestrow = 0; bestcol = 0;
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

	return BoardPosition {(int)bestrow, (int)bestcol};
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

	if (getNumHidden(game) < 20 && analytics.possibilities.logTotalCombinations < std::log(100000)) {
		//Use brute force to get perfect play
		BoardPosition bestMove = runBruteForce(game);

		Move out;
		out.flag = false;
		out.row = bestMove.row;
		out.col = bestMove.col;
		return out;
	}

	BoardPosition bestMove = runGS(game);

	Move out;
	out.flag = false;

	out.row = bestMove.row;
	out.col = bestMove.col;

	return out;
}

void Solver::clear() {
	queue.clear();
}