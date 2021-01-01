#pragma once
#include <iostream>
#include <vector>
#include <memory>
#include "Position.h"

			/********************************************/
			/************Part of task A.c****************/
			/***************Search Maze******************/
			/********************************************/

class Solution
{
public:
	Solution();
	//~Solution() {};
	std::shared_ptr<std::vector<std::unique_ptr<Position>>> getSolution();
	//std::vector<std::unique_ptr<Position>> getSolution1();
	void setSolution(std::unique_ptr<Position> s);
	int SolutionSize();


protected:
	std::vector<std::unique_ptr<Position>> m_solution;
	std::shared_ptr<std::vector<std::unique_ptr<Position>>> mShearedSolve;
	//Explain to myself about sheardPointer here
	//smart point is there in order to return it to function. 
	//It can not return as a vector value because it has a value of a pointer in it  ...
};

