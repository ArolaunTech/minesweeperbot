# Minesweeper bot

I don't have a name for this currently (though maybe I will eventually). Right now it's fairly strong, but not strong enough.

## Current win rate

91.26% &plusmn; 0.18% (100 000 trials) on Beginner **(9x9, 10 mines)**

78.10% &plusmn; 0.81% (10 000 trials) on Intermediate **(16x16, 40 mines)**

38.54% &plusmn; 0.95% (10 000 trials) on Expert **(16x30, 99 mines)**

## How it works

The solver uses a number of techniques to solve Minesweeper boards.

### Trivially safe / unsafe cells

In order to avoid doing expensive probability calculations every move, the solver checks how many flags and unrevealed cells are around every number on the board.

For example, the `5` in the image below had 5 empty squares and no flags around it. To satisfy the `5`, we need to flag every cell around it:

![All Mines](https://github.com/ArolaunTech/minesweeperbot/blob/main/imgs/5-flagged.png)

This then satisfies the `2` below it, allowing us to reveal a `3`:

![Satisfied](https://github.com/ArolaunTech/minesweeperbot/blob/main/imgs/2-rev.png)

This will get us pretty far especially on lower-difficulty boards.

### Mine Probabilities

Once we run out of trivially safe/unsafe cells to click/flag, we calculate the probability that each cell on the board is a mine by going through every possible combination of mines. We then click the safest cell.

## Build instructions

Run `sh build.sh` to compile the program. CMake is required.