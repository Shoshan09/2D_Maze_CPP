#include "Solution.h"


Solution::Solution() 
{	//Initializing the pointer to the array of the solution
	mShearedSolve = std::shared_ptr<std::vector<std::unique_ptr<Position>>>(&m_solution);
}
std::shared_ptr<std::vector<std::unique_ptr<Position>>> Solution::getSolution() 
{ 
	return mShearedSolve; 
}

void Solution::setSolution(std::unique_ptr<Position> s) {
	m_solution.push_back(std::move(s)); //movee is used for move the unique_ptr
}

int Solution::SolutionSize() {
	return m_solution.size();
}

//std::vector<std::unique_ptr<Position>> Solution::getSolution1() {
//	return m_solution;
//}

