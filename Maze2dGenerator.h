#pragma once
#include<iostream> 
#include <vector>
#include<string>		/*  for string class  */ 
#include <chrono>		/*  libraty  for time and date	*/
#include <stdlib.h>     /*	for rand metoed	 */
#include <algorithm>    // std::swap
#include "Maze2d.h"

			/********************************************/
			/************Part of task A.a****************/
			/***************Create Maze******************/
			/********************************************/

/*
**A class that contains a maze object
**and produces walls for it by the inheriting classes
**
**A maze is created  with 5*5 if there is no  given w and h
*/

class Maze2dGenerator //Abstract class
{
public:
	Maze2dGenerator();
	~Maze2dGenerator(){}
	Maze2dGenerator* operator = (Maze2dGenerator* mm); //as copy constractor


	virtual Maze2d& generate(int width, int height) =0 ; //Abstract method 
	virtual std::string measureAlgorithmTime(int width, int height);

	inline Maze2d getMaze() { return m; }

	//Position Can not exist in a class that is not a maze. Because the maze produces it 
	//and any way to get to it is only from the maze 
	//so it can not be returned as an independent type like: 
	//inline Position getPositions() { return m.Positions(); }

	//get
	inline int getStartP() { return m.Positions().getStartPosition(); }
	inline int getCurrP() { return m.Positions().getCuurPosition(); }
	inline int getEndP() { return m.Positions().getEndPosition(); }
	inline int getCurrPositionX() { return m.getCurrX(); }
	inline int getCurrPositionY() { return m.getCurrY(); }




	//maybe to delete
	//Maze2dGenerator(Position& p) { m.Positions() = p; }


	//set 
	inline void setStartP(int start){ m.Positions().StartPosition(start); }
	inline void setCurrP(int curr) { m.Positions().CuurPosition(curr); }
	inline void setEndP(int end) { m.Positions().EndPosition(end); }
	inline void setX(int x) { m.Positions().setX(x); }
	inline void setY(int y) { m.Positions().setY(y); }

protected:
	Maze2d m;
};



///***-------------------General maze----------------------------***///
class simpleMaze2dGenerator : public Maze2dGenerator
{
public:
	simpleMaze2dGenerator();
	//simpleMaze2dGenerator(const simpleMaze2dGenerator& m);
	virtual Maze2d& generate(int width, int height);

};


///***-------------------My maze by DFS----------------------------***///
/*DFS CODE In general, according to which the algorithm was built:
---------------------------------------------------------------------
map = [][] # with full of walls
for each (i, j) cell in map:
  DFS(i,j)

def DFS(i, j):
  if i, j is out of boundary:
	return
  if (i,j) is visited:
	return
  if number of visited neighbors of (i, j) > 1:
	return
  mark cell as visited
  for each neighbor (ni, nj) of cell(random order):
	DFS(ni, nj)
*/class MyMaze2dGenerator : public Maze2dGenerator
{
public:
	MyMaze2dGenerator();
	virtual Maze2d& generate(int width, int height);
	void generate_maze(int i, int j);
	int visitedAmount(int i, int j);
	void rendomWalls(int a[], int size);

private:
	//4 options to visit the cube :
	//  0-right , 1-left , 2-up, 3-down. 
	int visitOrderCell[4] = { 0,1,2,3 };
};


