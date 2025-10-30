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

bool isUnrelevantRadius(Game& game, int row, int col, int radius) {
	if (!game.cellIsState(row, col, CELL_HIDDEN)) return false;

	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	bool neighborless = true;
	for (int dr = -radius; dr <= radius; dr++) {
		for (int dc = -radius; dc <= radius; dc++) {
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

bool isUnrelevant(Game& game, int row, int col) {
	return isUnrelevantRadius(game, row, col, 1);
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

SearchResult Solver::search(
	Game& game, 
	std::vector<int> allowedIndices,
	std::vector<bool> allowedClicks,
	int totalCombinations, 
	std::vector<BoardPosition>& hiddenCells, 
	std::vector<std::vector<std::size_t> >& hiddenCellIndex,
	std::vector<std::vector<bool> >& mineCombinations,
	std::unordered_map<std::pair<std::vector<int>, std::vector<bool> >, SearchResult, TranspositionHash>& transpositionTable
) {
	bruteForceCalls++;

	if (transpositionTable.contains(std::pair<std::vector<int>, std::vector<bool> >(allowedIndices, allowedClicks))) {
		return transpositionTable[std::pair<std::vector<int>, std::vector<bool> >(allowedIndices, allowedClicks)];
	}

	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();

	std::size_t numHiddenCells = hiddenCells.size();

	std::size_t bestClick = 0;
	int bestWins = 0;
	bool forcedClick = false;

	std::vector<bool> dead(numHiddenCells, true);
	std::vector<int> potentialwinsarr(numHiddenCells);
	std::vector<int> clickorder(numHiddenCells);
	bool alldead = true;

	for (std::size_t click = 0; click < numHiddenCells; click++)
		clickorder[click] = click;

	for (std::size_t click = 0; click < numHiddenCells; click++) {
		if (!allowedClicks[click]) continue;

		BoardPosition clickPosition = hiddenCells[click];

		bool mine = false;
		int cellnumpossibilities = 0;
		int potentialwins = 0;

		std::vector<bool> seen(9);

		for (int i : allowedIndices) {
			if (mineCombinations[i][click]) {
				mine = true;
				continue;
			}

			potentialwins++;

			int cellnumber = 0;

			for (int dr = -1; dr <= 1; dr++) {
				for (int dc = -1; dc <= 1; dc++) {
					if (dr == 0 && dc == 0) continue;

					int newrow = clickPosition.row + dr;
					int newcol = clickPosition.col + dc;

					if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) continue;

					CellState other = game.getCell(newrow, newcol);

					if (other == CELL_FLAG) cellnumber++;

					if (other == CELL_HIDDEN && mineCombinations[i][hiddenCellIndex[newrow][newcol]]) {
						cellnumber++;
					}
				}
			}

			if (!seen[cellnumber]) {
				seen[cellnumber] = true;
				cellnumpossibilities++;
			}
		}

		potentialwinsarr[click] = potentialwins;

		if (cellnumpossibilities > 1) {
			dead[click] = false;
			alldead = false;
		}

		if (!mine) {
			bestClick = click;
			forcedClick = true;
			alldead = false;
			dead[click] = false;
			break;
		}
	}

	if (alldead) {
		SearchResult out;
		out.wins = 1;
		out.bestClick = 0;

		return out;
	}

	std::sort(clickorder.begin(), clickorder.end(), [potentialwinsarr](int a, int b) {
		return potentialwinsarr[a] > potentialwinsarr[b];
	});

	for (const int& click : clickorder) {
		if (!allowedClicks[click]) continue;
		if (forcedClick && click != bestClick) continue;
		if (dead[click]) continue;

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

		if (potentialwins <= bestWins) break;

		int wins = 0;
		for (int i = 0; i < 9; i++) {
			std::size_t numRemainingPossibilities = indices[i].size();

			if (numRemainingPossibilities == 0) continue;
			if (numRemainingPossibilities == 1) {
				wins++;
				continue;
			}

			int searchwins = search(
				game, 
				indices[i], 
				newAllowedClicks, 
				totalCombinations, 
				hiddenCells, 
				hiddenCellIndex, 
				mineCombinations,
				transpositionTable
			).wins;

			potentialwins += searchwins - numRemainingPossibilities;

			if (potentialwins <= bestWins) break;

			wins += searchwins;
		}

		if (wins > bestWins) {
			bestClick = click;
			bestWins = wins;
		}
	}

	SearchResult out;
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

	std::unordered_map<std::pair<std::vector<int>, std::vector<bool> >, SearchResult, TranspositionHash> transpositionTable;

	SearchResult result = search(
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
	//IDEAS: 
	// - some sort of RAVE
	// - EMCTS
	// - whatever I read
	// - limited tree search with ~1000 mine configs
	std::size_t numrows = game.getRows();
	std::size_t numcols = game.getCols();
	int totalmines = game.getTotalMines();

	AnalyzeResult analytics = analyze(game);

	std::vector<std::vector<double> > probabilities = analytics.probabilities;
	PossibilitiesResult possibilities = analytics.possibilities;

	int numflags = getNumFlags(game);

	std::size_t numRelevantCellGroups = possibilities.relevantCells.size();

	std::size_t bestrow = 0, bestcol = 0;

	/*
	//Weird RAVE test
	std::vector<BoardPosition> moves;
	for (std::size_t i = 0; i < numrows; i++) {
		for (std::size_t j = 0; j < numcols; j++) {
			if (!game.cellIsState(i, j, CELL_HIDDEN)) continue;
			if (isUnrelevantRadius(game, i, j, 2)) continue;

			moves.push_back(BoardPosition {(int)i, (int)j});
		}
	}

	std::size_t nummoves = moves.size();

	std::vector<double> score(nummoves);
	std::vector<int> plays(nummoves);
	int totalplays = 0;

	double initconfidence = 10;
	for (std::size_t i = 0; i < nummoves; i++) {
		score[i] = initconfidence * (1 - probabilities[moves[i].row][moves[i].col]);
		plays[i] = initconfidence;
	}

	for (int i = 0; i < 100 * nummoves; i++) {
		double selector = randfloat(0, 1);

		int possIndex = 0;
		while (selector >= 0) {
			if (!possibilities.valid[possIndex]) {
				possIndex++;
				continue;
			}

			selector -= std::exp(possibilities.logCombinations[possIndex] - possibilities.logTotalCombinations);
			possIndex++;
		}

		possIndex--;

		std::vector<std::size_t> possibility = possibilities.possibilities[possIndex];
		std::vector<std::vector<bool> > mines(numrows, std::vector<bool>(numcols));

		for (std::size_t row = 0; row < numrows; row++) {
			for (std::size_t col = 0; col < numcols; col++) {
				mines[row][col] = game.cellIsState(row, col, CELL_FLAG);
			}
		}

		int irrelevantMines = totalmines - numflags;

		for (std::size_t group = 0; group < numRelevantCellGroups; group++) {
			std::size_t numcells = possibilities.relevantCells[group].size();
			std::size_t nummines = possibility[group];

			irrelevantMines -= nummines;

			for (std::size_t mine = 0; mine < nummines; mine++) {
				int randrow, randcol;

				do {
					int randidx = randint(0, numcells - 1);

					randrow = possibilities.relevantCells[group][randidx].row;
					randcol = possibilities.relevantCells[group][randidx].col;
				} while (mines[randrow][randcol]);

				mines[randrow][randcol] = true;
			}
		}
		
		for (int mine = 0; mine < irrelevantMines; mine++) {
			int randrow, randcol;

			do {
				randrow = randint(0, numrows - 1);
				randcol = randint(0, numcols - 1);
			} while (!isUnrelevant(game, randrow, randcol) || mines[randrow][randcol]);

			mines[randrow][randcol] = true;
		}

		for (int j = 0; j < 1; j++) {
			std::vector<std::vector<bool> > revealed = game.getRevealed();
			std::vector<int> movesMade;

			bool loss = false;
			for (int move = 0; move < 50; move++) {
				int bestrandidx = 0;
				double bestrandscore = 0;

				for (int movetry = 0; movetry < nummoves; movetry++) {
					int randidx = movetry;

					if (revealed[moves[randidx].row][moves[randidx].col]) continue;

					double currscore = score[randidx] / plays[randidx] + 1.3 * std::sqrt(std::log(totalplays) / plays[randidx]);

					if (currscore > bestrandscore) {
						bestrandidx = randidx;
						bestrandscore = currscore;
					}
				}

				BoardPosition position = moves[bestrandidx];

				if (mines[position.row][position.col]) {
					loss = true;
				}

				if (!revealed[position.row][position.col]) movesMade.push_back(bestrandidx);

				revealed[position.row][position.col] = true;

				//if (randchance(1.0 / 2.0)) break;
				break;
			}

			for (std::size_t move = 0; move < movesMade.size(); move++) {
				plays[movesMade[move]]++;
				totalplays++;
			}

			if (loss) continue;

			std::vector<std::vector<bool> > flagged = game.getFlagged();

			bool moved;
			do {
				moved = false;

				for (std::size_t row = 0; row < numrows; row++) {
					for (std::size_t col = 0; col < numcols; col++) {
						if (!revealed[row][col]) continue;

						bool unflaggedmine = false;
						bool unrevealedsafe = false;

						for (int dr = -1; dr <= 1; dr++) {
							for (int dc = -1; dc <= 1; dc++) {
								if (dr == 0 && dc == 0) continue;

								int newrow = row + dr;
								int newcol = col + dc;

								if (newrow < 0 || newrow >= numrows) continue;
								if (newcol < 0 || newcol >= numcols) continue;

								if (mines[newrow][newcol]) {
									if (!flagged[newrow][newcol]) unflaggedmine = true;
								} else {
									if (!revealed[newrow][newcol]) unrevealedsafe = true;
								}
							}
						}

						//!unflaggedmine -> reveal all cells around this one
						//unflaggedmine && !unrevealedsafe -> flag all cells around this one

						if (unflaggedmine == unrevealedsafe) continue;

						for (int dr = -1; dr <= 1; dr++) {
							for (int dc = -1; dc <= 1; dc++) {
								if (dr == 0 && dc == 0) continue;

								int newrow = row + dr;
								int newcol = col + dc;

								if (newrow < 0 || newrow >= numrows) continue;
								if (newcol < 0 || newcol >= numcols) continue;

								if (revealed[newrow][newcol] || flagged[newrow][newcol]) continue;

								moved = true;

								if (unflaggedmine) {
									flagged[newrow][newcol] = true;
								} else {
									revealed[newrow][newcol] = true;
								}
							}
						}
					}
				}
			} while (moved);

			std::vector<double> weightscells = {
				-0.0028239844516102157,
				-0.0002665864378393504, 
				0.010053411630445168, 
				0.00787855859419873, 
				0.0021288558858229673, 
				0.025608179166585483, 
				-0.008110171201601629, 
				-0.0075273725447626135, 
				-0.009202722549098787, 
				0.0034025193071922696, 
				0.0005806418980248698, 
				0.00803225396542895
			};
			std::vector<double> weightspoly = {
				0.006852606748328438,
				0.002992242857208602, 
				0.0003433411838952726, 
				0.01899811373350083, 
				0.020957357631834643, 
				0.0053231122830747005
			};

			double finalsum = 0;
			for (std::size_t row = 0; row < numrows; row++) {
				for (std::size_t col = 0; col < numcols; col++) {
					finalsum += weightscells[game.getCell(row, col)];
				}
			}

			finalsum *= 480.0 / numrows / numcols;

			double finalscore = 0;
			for (std::size_t term = 0; term < weightspoly.size(); term++) {
				finalscore += std::pow(finalsum, term) * weightspoly[term];
			}

			finalscore = 1 - 1 / (1 + std::exp(finalscore));

			//if (finalscore < .1) finalscore *= 5;
			//else if (finalscore > .9) finalscore = 1 - (1 - finalscore) * 5;
			//else finalscore = .5;

			for (std::size_t move = 0; move < movesMade.size(); move++) {
				score[movesMade[move]] += finalscore;
			}
		}
	}

	//logVector(plays);
	//std::cout << "plays\n";
	//logVector(score);
	//std::cout << "score\n";

	std::size_t bestidx = 0;
	for (std::size_t i = 0; i < nummoves; i++) {
		if (plays[i] > plays[bestidx]) {
			bestidx = i;
		}
	}

	bestrow = moves[bestidx].row;
	bestcol = moves[bestidx].col;
	*/

	int iters = 1000;
	int n = 100;

	std::vector<GameTree> population(n);

	for (int i = 0; i < iters; i++) {
		double selector = randfloat(0, 1);

		int possIndex = 0;
		while (selector >= 0) {
			if (!possibilities.valid[possIndex]) {
				possIndex++;
				continue;
			}

			selector -= std::exp(possibilities.logCombinations[possIndex] - possibilities.logTotalCombinations);
			possIndex++;
		}

		possIndex--;

		std::vector<std::size_t> possibility = possibilities.possibilities[possIndex];
		std::vector<std::vector<bool> > mines(numrows, std::vector<bool>(numcols));

		for (std::size_t row = 0; row < numrows; row++) {
			for (std::size_t col = 0; col < numcols; col++) {
				mines[row][col] = game.cellIsState(row, col, CELL_FLAG);
			}
		}

		int irrelevantMines = totalmines - numflags;

		for (std::size_t group = 0; group < numRelevantCellGroups; group++) {
			std::size_t numcells = possibilities.relevantCells[group].size();
			std::size_t nummines = possibility[group];

			irrelevantMines -= nummines;

			for (std::size_t mine = 0; mine < nummines; mine++) {
				int randrow, randcol;

				do {
					int randidx = randint(0, numcells - 1);

					randrow = possibilities.relevantCells[group][randidx].row;
					randcol = possibilities.relevantCells[group][randidx].col;
				} while (mines[randrow][randcol]);

				mines[randrow][randcol] = true;
			}
		}
		
		for (int mine = 0; mine < irrelevantMines; mine++) {
			int randrow, randcol;

			do {
				randrow = randint(0, numrows - 1);
				randcol = randint(0, numcols - 1);
			} while (!isUnrelevant(game, randrow, randcol) || mines[randrow][randcol]);

			mines[randrow][randcol] = true;
		}

		std::vector<std::vector<int> > cellnumbers(numrows, std::vector<int>(numcols));
		for (std::size_t row = 0; row < numrows; row++) {
			for (std::size_t col = 0; col < numcols; col++) {
				for (int dr = -1; dr <= 1; dr++) {
					for (int dc = -1; dc <= 1; dc++) {
						if (dr == 0 && dc == 0) continue;

						int newrow = row + dr;
						int newcol = col + dc;

						if (newrow < 0 || newrow >= numrows || newcol < 0 || newcol >= numcols) continue;

						if (mines[newrow][newcol]) cellnumbers[row][col]++;
					}
				}
			}
		}

		int idx1 = randint(0, n - 1);
		int idx2 = 0;

		do {
			idx2 = randint(0, n - 1);
		} while (idx1 == idx2);

		GameTree tree1 = population[idx1];
		GameTree tree2 = population[idx2];

		double score1 = 0, score2 = 0;
		bool lost1 = false, lost2 = false;

		std::vector<std::vector<bool> > revealed = game.getRevealed();

		while (true) {
			BoardPosition move = tree1.move;

			if (mines[move.row][move.col]) {
				lost1 = true;
				break;
			}

			revealed[move.row][move.col] = true;

			if (tree1.children.empty()) {
				//Evaluate
				score1 = 1;
				break;
			}

			tree1 = tree1.children[cellnumbers[move.row][move.col]];
		}

		revealed = game.getRevealed();

		while (true) {
			BoardPosition move = tree2.move;

			if (mines[move.row][move.col]) {
				lost2 = true;
				break;
			}

			revealed[move.row][move.col] = true;

			if (tree2.children.empty()) {
				//Evaluate
				score2 = 1;
				break;
			}

			tree2 = tree2.children[cellnumbers[move.row][move.col]];
		}

		if (lost1) score1 = -1;
		if (lost2) score2 = -1;

		if (score1 == score2) continue;
		if (score1 > score2) {

		} else {

		}
	}

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

	if (getNumHidden(game) < 22 && analytics.possibilities.logTotalCombinations < std::log(100000)) {
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