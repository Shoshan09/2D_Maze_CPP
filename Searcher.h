#pragma once
#include "Solution.h"
#include "Searchable.h"

			/********************************************/
			/************Part of task A.c****************/
			/***************Search Maze******************/
			/********************************************/

class Searcher //Pure Abstract class 
{
public:
	virtual Solution search(Searchable& s) = 0;
	virtual int getNumberOfNodeMaze2d() = 0;  //how much has allrady invaloated 
};




