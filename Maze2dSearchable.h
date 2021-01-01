#pragma once
#include<iostream> 
#include<string>		/*  for string class  */ 
#include "Position.h"
#include "Maze2dGenerator.h"
#include "Searchable.h"

			/********************************************/
			/************Part of task A.c****************/
			/***************Search Maze******************/
			/********************************************/


/* Algoritem Search For Maze:
** I chose to use a maze2d object because I check what the situation is in the maze. 
** in the main will select the type from inhariting class and initialize to the size and the size 
** the maze thet will be returned is the maze itself and will be dropped into this class Maze2dSearchable

**Could also be realized with a constructor containing class Maze2dGenerator, I did not choose to do so.
*/

class Maze2dSearchable : public Searchable
{
public:
	Maze2dSearchable(Maze2d& maze);
	Maze2dSearchable* operator = (Position* p);

	////implamantion of Searchable mathodes
	virtual std::vector<std::string> getpossiblemoves(Position& p) override;
	virtual Position getpositions() override;
	virtual void setNewPosition(int x, int y) override;
	virtual std::string getClass() { return "Maze2dSearchable";}
	Maze2d& getMaze() { return m; }


protected:
	Maze2d m;
	Position p;			/*To initialize the user's points from the maze --us and chang it as searching for a path */
};




//class Maze2dSearchable : public Searchable
//{
//public:
//	Maze2dSearchable(Maze2d& maze);
//	Maze2dSearchable* operator = (Position* p);
//
//	////implamantion of Searchable mathodes
//	virtual std::vector<std::string> getpossiblemoves(Position& p) override;
//	virtual Position getpositions() override;
//	virtual void setNewPosition(int x, int y) override;
//	//virtual T& getGame() { return m; }
//
//protected:
//	Maze2d m;
//	Position p;			/*To initialize the user's points from the maze --us and chang it as searching for a path */
//};

