#include "Demo.h"
#include "CommonSearcher.h"


void Demo::run() {
	//a:
	MyMaze2dGenerator m4; 
	Maze2d& maze = m4.generate(30, 30);
	//b:
	std::cout << maze << std::endl;
	//c:
	Maze2dSearchable e(maze);
	BFS bfs;
	Solution b = bfs.search(e);
	//d:
	//AStar astar;
	//Solution as = astar.search(e);
	//astar.setSolve(e, b);			//set in the maze
	//astar.printAStarInMaze(e);
	//e:
	//std::cout << " Number Of Node Maze2d BFS: " << bfs.getNumberOfNodeMaze2d() << std::endl;
	//std::cout << " Number Of Node Maze2d Astar: " << astar.getNumberOfNodeMaze2d() << std::endl
}
