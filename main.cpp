#include <iostream>
#include <utility>
using namespace std;
#include "Maze2dGenerator.h"
#include "TestMazeGenerator.h"
#include "Maze2dSearchable.h"
#include "CommonSearcher.h"
#include "Demo.h"


#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3


int main(void)
{
	
	TestMazeGenerator t;
	simpleMaze2dGenerator m2;
	t.testMazeGenerator(m2);
	//MyMaze2dGenerator m3;
	//t.testMazeGenerator(m3);

	//main for searches:
	MyMaze2dGenerator m4; //create by defulte 5*5
	Maze2d& maze = m4.generate(15, 15);

	Maze2dSearchable e(maze);
	auto p = e.getpositions();
	std::vector<std::string> s = e.getpossiblemoves(p);
	BFS bfs;
	bool a[2];
	//bfs.isValid(e, a, 0, 0);
	Solution b = bfs.search(e);
	//Solution s = bfs.search(e);
	AStar astar;
	astar.aStarSearch(e, int grid[][COL], Pair src, Pair dest);

	// Destination is the left-most top-most corner 
	

	return 0;
}