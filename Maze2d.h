#pragma once
#include<iostream> 
#include <vector>
#include <memory>			// for samrt pointers
#include<string>			/*  for string class  */ 
#include <fstream>			 //ifstream
#include <map>			
#include "Position.h"


			/********************************************/
			/************Part of task A.a****************/
			/***************Create Maze******************/
			/********************************************/



/*
** A class that produces a matrix,
** 0 is an an empty space
** 1 is a wall.
** 2 is for the passage visited (for DFS)
** 9 is for a starting position
** 10 is for a goal posithion
**
** The walls are realized In a class that contains a maze-type object
*/

class Maze2d
{
public:
	Maze2d(int width, int height);
	Maze2d();
	Maze2d* operator = (Maze2d* m); //as copy constractor
	Maze2d& operator = (Maze2d& m) {
		_maze = m._maze;
		sidesMazeIndexs = m.sidesMazeIndexs;
		_tempArray = m._tempArray;
		w = m.w;
		h = m.h;
		_mazPositions = m._mazPositions;
		startIndex = m.startIndex;
		currIndex = m.currIndex;
		endIndex = m.endIndex;
		return *this;
	}

	//get
	inline auto getMaze(int w, int h) { return _maze[w][h]; }
	inline int getEnd() { return  endIndex; }
	inline int getH() { return h; }
	inline int getW() { return w; }
	inline int getStart() { return startIndex; }	
	int getCurrMazeIndex();
	inline int getCurrX() { return currX; }
	inline int getCurrY() { return currY; }
	inline auto getSizeMaze() { return _maze.size(); }
	inline auto rowSize() { return _maze[0].size(); }
	inline Position Positions() { return _mazPositions; }
	std::vector < std::vector<int>> getData() { return _maze; }


	inline int getIndex(int r,int c) { 
		for (int i = 0; i < w; i++)
			for (int j = 0; j < h; j++)
				if (i == r && j == c)
					return _maze[i][j];
	}


	//functions as asked in the project
	std::string getGoalPosition();
	std::vector<std::string> getPossibleMoves(Position& p);
	Position getStartPosition();

	//set
	void setMaze(int i, int j, int num);
	void setCurrPosition(int curr);

	//dont need it apleyed in operator <<
	void drowMaze();

	//Realization of printing
	friend std::ostream& operator << (std::ostream& os, const Maze2d& p) {
		char Wall = char(219);
		os << "S is for Starting Point\nG is for Goal Point\nThe Maze:  \n\n";
	
		 
		for (int j = 0; j <= p.w; j++) 
			os << Wall;
		os << Wall;
		os << '\n';

		for (int i = 0; i < p.w; i++) {
			os << Wall;
			for (int j = 0; j < p.h; j++) {
				if (p._maze[i][j] == 0)
					os << " ";
				if (p._maze[i][j] == 1)
					os << Wall;
				if (p._maze[i][j] == 9)
					os << "S";
				if (p._maze[i][j] == 10)
					os << "G";
				if (p._maze[i][j] == 2) {  //path 
					os << "*";
				}
			}
			os << Wall;
			os << '\n';
		}

		for (int j = 0; j <= p.w; j++)
			os << Wall;
		os << Wall;
		os << '\n';
		return os;
	}



protected:
	int w, h;										/* width , height*/
	std::vector<std::vector<int>> _maze;			/* Stores maze	*/
	std::vector<int> sidesMazeIndexs;				/* Stores rows of maze exits*/
	std::vector<std::vector<int>>  _tempArray;		/* Stores rows of maze to be allocated into maze 2-dimensional vector */
	Position _mazPositions;							/*To initialize the user's points */
	double _cost;									/* Number of visit for shearching */
	int startIndex, currIndex, endIndex;			/* Stores cube number of maze positions		*/
	int currX, currY;
	char buffer;									/* To store char elements from external text-file	 */
};
