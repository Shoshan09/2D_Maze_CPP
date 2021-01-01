#pragma once
#include "Searcher.h"
#include "Maze2dSearchable.h"
#include <functional>
#include <queue>
#include <iostream>
#include <utility>


#define ROW 9
#define COL 10

using namespace std;

			/********************************************/
			/************Part of task A.c****************/
			/***************Search Maze******************/
			/********************************************/


class CommonSearcher : public Searcher
{
public:
	CommonSearcher();

	////implamantion of Searcher mathodes
	virtual Solution search(Searchable& s) = 0 ; //abstract mathode
	virtual int getNumberOfNodeMaze2d() { return m_evaluateNodes; }

	const Position& popOpenList() { m_evaluateNodes++; return m_openList.top(); }
	void print_queue();

protected:
	int m_evaluateNodes;
	std::priority_queue<Position> m_openList;
	int row[4], col[4];  
	// m_evaluateNodes is Distance from start to finish

};


//
/////***------------------BFS----------------------------***///
///*Pseudocode BFS:
//
//BFS (G, s)											//Where G is the graph and s is the source node
//	  let Q be queue.
//	  Q.enqueue( s )								//Inserting s in queue until all its neighbour vertices are marked.
//
//	  mark s as visited.
//	  while ( Q is not empty)
//		 v  =  Q.dequeue( )							 //Removing that vertex from queue,whose neighbour will be visited now
//		  
//													 //processing all the neighbours of v
//		  for all neighbours w of v in Graph G
//			   if w is not visited
//						Q.enqueue( w )				 //Stores w in Q to further visit its neighbour
//						mark w as visited.
//*/
//


class BFS : public CommonSearcher
{
public:
	BFS() : CommonSearcher(){}
	Solution search(Searchable& s);
	Solution search(Maze2dSearchable& s);

	bool isValid(Maze2dSearchable& s, std::vector<std::vector<int>> visited, int row, int col);
	bool isValidEnd(Maze2dSearchable& s, std::vector<std::vector<int>> visited, int row, int col);
	void drowFormat(Maze2dSearchable& s, std::vector<std::vector<int>>  visited);
	//int bfsMaze(Maze2dSearchable& s, std::vector<std::vector<int>> current, int row, int col);

private: 
	std::vector<std::unique_ptr<Position>> prev;
	bool GoalFound;
};



///***------------------AStar----------------------------***///

class Cell: public CommonSearcher
{
public:
	Cell() : CommonSearcher() {};
	// Row and Column index of its parent 
	// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1 
	int parent_i, parent_j;
	// f = g + h  
	double f, g, h;
	// Creating a shortcut for int, int pair type 
	typedef std::pair <int, int> Pair;

	// Creating a shortcut for pair<int, pair<int, int>> type 
	typedef std::pair <double, pair<int, int>> pPair;
	Pair src = make_pair(8, 0);

	// Destination is the left-most top-most corner 
	Pair dest = make_pair(0, 0);
};

class AStar : public CommonSearcher
{

public:
	
	AStar() : CommonSearcher() {}
	class Cell 
	{
	public:
		Cell();
		// Row and Column index of its parent 
		// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1 
		int parent_i, parent_j;
		// f = g + h  
		double f, g, h;
		// Creating a shortcut for int, int pair type 
		typedef std::pair <int, int> Pair;

		// Creating a shortcut for pair<int, pair<int, int>> type 
		typedef std::pair <double, pair<int, int>> pPair;
		Pair src = make_pair(8, 0);

		// Destination is the left-most top-most corner 
		Pair dest = make_pair(0, 0);
	};
	
	

	static bool isValid(int row, int col);
	bool isUnBlocked(Maze2dSearchable& s, int row, int col);
	bool isDestination(int row, int col, Cell::Pair dest);
	double calculateHValue(int row, int col, Cell::Pair dest);
	void tracePath(Cell cellDetails[][COL], Cell::Pair dest);
	void aStarSearch(Cell cellDetails[][COL],Maze2dSearchable& s, Cell::Pair src, Cell::Pair dest);
	

	
};

