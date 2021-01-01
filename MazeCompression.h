#pragma once
#include "Maze2d.h"
#include "Solution.h"

				/********************************************/
				/************Part of task B.a****************/
				/***************Use Files *******************/
				/********************************************/

class MazeCompression
{
	MazeCompression(Maze2d& m);
	void save(const char* FileName);
	void saveSolutionToStock();
	Solution getExistSolution();

private:
	std::vector<std::vector<int>> my_maze;			/* Stores maze	*/
	int w, h;
};

