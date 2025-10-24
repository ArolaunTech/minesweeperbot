#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "board.h"
#include "game.h"
#include "solver.h"

int main() {
	int numrows = 9;
	int numcols = 9;
	int totalmines = 10;

	Solver solver;

	int wins = 0;
	int losses = 0;

	for (int games = 0; games < 100000; games++) {
		Game game(numrows, numcols, totalmines);

		solver.clear();
		while (true) {
			//std::cout << game.toString();

			//double error = 1.96 * std::sqrt((double)wins * losses / (wins + losses) / (wins + losses) / games);
			//error = 0.01 * std::round(10000.0 * error);

			//std::cout << wins << " " << losses << " " << 0.01 * std::round(10000.0 * wins / (wins + losses)) << "% +- " << error << "%\n";

			GameState currstate = game.getGameState();

			if (currstate == GAME_LOSS) {
				//std::cout << "Loss :(\n" /*<< game.toString()*/;
				losses++;
				break;
			}

			if (currstate == GAME_WIN) {
				//std::cout << "Win!\n";
				wins++;
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

		double error = 1.96 * std::sqrt((double)wins * losses / (wins + losses) / (wins + losses) / games);
		error = 0.01 * std::round(10000.0 * error);

		std::cout << wins << " " << losses << " " << 0.01 * std::round(10000.0 * wins / (wins + losses)) << "% +- " << error << "%\n";
	}

	std::cout << wins << " " << losses << " " << 0.01 * std::round(10000.0 * wins / (wins + losses)) << "%\n";
}