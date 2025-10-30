# Minesweeper bot

I don't have a name for this currently (though maybe I will eventually). Right now it's fairly strong, but not strong enough.

## Best win rate

91.40% &plusmn; 0.17% (100 000 trials) on Beginner **(9x9, 10 mines)**

77.29% &plusmn; 0.26% (100 000 trials) on Intermediate **(16x16, 40 mines)**

39.33% &plusmn; 0.30% (100 000 trials) on Expert **(16x30, 99 mines)**

## How it works

The solver uses a number of techniques to solve Minesweeper boards. I list them here in order of least to most expensive.

### Trivially safe / unsafe cells

In order to avoid doing expensive probability calculations every move, the solver checks how many flags and unrevealed cells are around every number on the board.

For example, the `5` in the image below had 5 empty squares and no flags around it. To satisfy the `5`, we need to flag every cell around it:

![All Mines](https://github.com/ArolaunTech/minesweeperbot/blob/main/imgs/5-flagged.png)

This then satisfies the `2` below it, allowing us to reveal a `3`:

![Satisfied](https://github.com/ArolaunTech/minesweeperbot/blob/main/imgs/2-rev.png)

This will get us pretty far especially on lower-difficulty boards.

### Mine Probabilities

Once we run out of trivially safe/unsafe cells to click/flag, we calculate the probability that each cell on the board is a mine by going through every possible combination of mines. This allows us to solve patterns like the `2 1`, `1 2 1`, or `1 2 2 1`. This usually reveals many more safe and unsafe cells which we can then click/flag. If we are in the early-mid game, we click the safest cell, but if we are in the late-game, we use...

### Brute force search

When the number of possible mine configurations are low (< 100 000), the solver uses brute force to calculate the perfect move from endgame positions. This adds about 0.7% to the winrate on Expert mode.

## Weaknesses

The solver currently struggles with early-game play if it does not get a 0 on the first move.

## Build instructions

Run `sh build.sh` to compile the program. CMake is required.