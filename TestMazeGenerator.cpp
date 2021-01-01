#include "TestMazeGenerator.h"


void TestMazeGenerator::testMazeGenerator(Maze2dGenerator& mg) {
	cout << mg.measureAlgorithmTime(15,15) << endl;
	 Maze2d& maze = mg.generate(15, 15);
	Position p = maze.getStartPosition();
	cout << p << endl;
	std::vector<std::string> moves = maze.getPossibleMoves(p);
	for (string move : moves) 
		cout << move << endl;
	cout << maze.getGoalPosition() << endl;
	cout << maze << endl;
}
