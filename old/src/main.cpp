#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>

#include "board.h"
#include "game.h"
#include "solver.h"
#include "debug.h"
#include "random.h"

int main() {
	int numrows = 16;
	int numcols = 30;
	int totalmines = 99;

	Solver solver;
	solver.bruteForceCalls = 0;

	int wins = 0;
	int losses = 0;

	std::vector<int> revealedDist(numrows * numcols);

	//std::ofstream log("log.txt");

	for (int games = 0; games < 100000; games++) {
		Game game(numrows, numcols, totalmines);

		solver.clear();
		while (true) {
			//std::cout << game.toString();

			//double error = 1.96 * std::sqrt((double)wins * losses / (wins + losses) / (wins + losses) / games);
			//error = 0.01 * std::round(10000.0 * error);

			//std::cout << wins << " " << losses << " " << 0.01 * std::round(10000.0 * wins / (wins + losses)) << "% +- " << error << "%\n";

			//std::vector<int> freqs(12);
			//for (int i = 0; i < numrows; i++) {
			//	for (int j = 0; j < numcols; j++) {
			//		freqs[game.getCell(i, j)]++;
			//	}
			//}

			//for (int i = 0; i < 12; i++) {
			//	log << freqs[i] << " ";
			//}
			//log << "\n";

			//log << solver.analyze(game).possibilities.logTotalCombinations << " ";

			GameState currstate = game.getGameState();

			if (currstate == GAME_LOSS) {
				//std::cout << "Loss :(\n" /*<< game.toString()*/;
				losses++;

				int revealed = 0;
				for (int i = 0; i < numrows; i++) {
					for (int j = 0; j < numcols; j++) {
						CellState cell = game.getCell(i, j);
						if (cell != CELL_HIDDEN && cell != CELL_FLAG && cell != CELL_MINE) {
							revealed++;
						}
					}
				}

				revealedDist[revealed]++;

				//log << "loss\n";
				break;
			}

			if (currstate == GAME_WIN) {
				//std::cout << "Win!\n";
				wins++;

				//log << "win\n";
				break;
			}

			/*std::vector<std::vector<double> > probabilities = solver.getMineProbabilities(game);
	
			for (int i = 0; i < numrows; i++) {
				for (int j = 0; j < numcols; j++) {
					int roundedProbability = (int)std::round(100 * probabilities[i][j]);
	
					if (roundedProbability < 10) {
						std::cout << "  " << roundedProbability << " ";
					} else if (roundedProbability == 100) {
						std::cout << "100 ";
					} else {
						std::cout << " " << roundedProbability << " ";
					}
				}
				std::cout << "\n";
			}*/

			Move move = solver.getBestMove(game);

			if (move.flag) game.flag(move.row, move.col); else game.click(move.row, move.col);
		}
		

		if /*(games % 100 == 99)*/(true) {
			double error = 1.96 * std::sqrt((double)wins * losses / (wins + losses) / (wins + losses) / games);
			error = 0.01 * std::round(10000.0 * error);

			std::cout << wins << " " << losses << " " << 0.01 * std::round(10000.0 * wins / (wins + losses)) << "% +- " << error << "%\n";
			std::cout << (double)solver.bruteForceCalls / (wins + losses) << "\n";
		}
	}

	//log.close();

	std::cout << wins << " " << losses << " " << 0.01 * std::round(10000.0 * wins / (wins + losses)) << "%\n";
	logVector(revealedDist);
	std::cout << "\n";
}