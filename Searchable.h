#pragma once
#include<iostream>
#include<string>		/*  for string class  */ 
#include "Maze2dGenerator.h"
#include "Position.h"

			/********************************************/
			/************Part of task A.c****************/
			/***************Search Maze******************/
			/********************************************/

class Searchable  //Pure Abstract class 
{
public:
	virtual std::vector<std::string> getpossiblemoves(Position& p) = 0 ;
	virtual Position getpositions() = 0 ;
	virtual void setNewPosition(int x, int y) = 0;
	virtual std::string getClass() = 0;

};

//
//template <class T>
//class Maze2dSearchable : public Searchable<Maze2d>
//{
//public:
//	Maze2dSearchable(Maze2d& maze);
//	Maze2dSearchable* operator = (Position* p);
//
//	////implamantion of Searchable mathodes
//	virtual std::vector<std::string> getpossiblemoves(Position& p) override;
//	virtual Position getpositions() override;
//	virtual void setNewPosition(int x, int y) override;
//	//virtual Maze2d& getGame() { return m; }
//	virtual T& getGame() { return m; }
//
//
//
//
//protected:
//	Maze2d m;
//	Position p;			/*To initialize the user's points from the maze --us and chang it as searching for a path */
//};



