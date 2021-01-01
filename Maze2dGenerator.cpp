#include "Maze2dGenerator.h"


Maze2dGenerator::Maze2dGenerator() {}

Maze2dGenerator* Maze2dGenerator::operator = (Maze2dGenerator* mm) {
	m  = mm->m; //send to opertor =  of class Maze2d 
	return this;
}


std::string Maze2dGenerator::measureAlgorithmTime(int width, int height) {
	auto start = std::chrono::high_resolution_clock::now();
	//measure time for function :
	generate(width, height);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> duration = end - start;

	float ms = duration.count() * 1000.0f;
	std::string time = "Time took " + std::to_string(ms) + "ms";

	return time;
}

//void Maze2dGenerator::drowFormat() {
//	int num = 1;
//	for (int i = 0; i < m.getWidth(); i++)
//	{
//		std::cout << "  ";
//		for (int j = 0; j < m.getHeight(); j++)
//		{
//			std::cout << " " << num++ << "  ";
//		}
//		std::cout << " ";
//		std::cout << '\n';
//	}
//}
//
//void Maze2dGenerator::drow() {
//	for (int i = 0; i < m.getWidth(); i++) {
//		std::cout << "  ";
//		for (int j = 0; j < m.getHeight(); j++)
//		{
//			if (m.getMaze(i,j) == ' ')
//				std::cout << " ";
//			else
//				std::cout << m.getMaze(i, j);
//			std::cout << "  ";
//		}
//		std::cout << '\n';
//	}
//}



///***-------------------General maze----------------------------***///
simpleMaze2dGenerator::simpleMaze2dGenerator() : Maze2dGenerator() {}
//need Fixing !!!  
Maze2d& simpleMaze2dGenerator::generate(int width, int height) {
	//Create a Maze in the Size the User requsted and set it as an object of the class!!
	Maze2d m11(width, height);
	m = m11;

	/*Terms for testing:
	**  Wall =1;
	**  Free = 0;
	**  PlayerPosition = 9;
	**  GoalPosition = 10;
	*/

	//Set inside Walls Remdomly 
	//Definition of Chacking : 
	//			 1. if its Free 
	//			 2. if its not bloking the PlayerPosition or the  GoalPosition
	//			 3. Hold a temporary array that contains the walls of the row above the curr row 
	//				and check that the one above it is also a wall 

	int rowRow;
	int colum;
	////Random selection point
	//colum = rand() % m.getW();
	//rowRow = rand() % m.getH();

	for (int i = 1; i < m.getW() / 2; i++)
		for (int ii = 1; ii < m.getH() / 2; ii++) {
			//Random selection point
			colum = rand() % m.getW();
			rowRow = rand() % m.getH();

			//Check if colum, rowRow is Within the matrix
			if (colum == 0 || rowRow == 0 || colum == m.getH() || rowRow == m.getW()) break;

			//Checking by Definition and Set Walls in the Maze
			if (m.getMaze(colum, rowRow) == 0 && (m.getMaze(colum - 1, rowRow) != 9 || m.getMaze(colum + 2, rowRow) != 9 || m.getMaze(colum - 1, rowRow + 1) != 9 || m.getMaze(colum, rowRow - 2) != 9))
				if (m.getMaze(colum - 1, rowRow) != 10 || m.getMaze(colum + 2, rowRow) != 10 || m.getMaze(colum - 1, rowRow + 1) != 10 || m.getMaze(colum, rowRow - 2) != 10)
					m.setMaze(colum, rowRow, 1);
		}
	return m;
}


///***-------------------my maze by DFS----------------------------***///
MyMaze2dGenerator::MyMaze2dGenerator() : Maze2dGenerator() {}

int MyMaze2dGenerator::visitedAmount(int i, int j) {
	int direct[][2] = { {1,0}, {-1,0}, {0,1}, {0,-1} };
	int count = 0;
	for (int k = 0; k < 4; ++k)
	{
		int ni = i + direct[k][0];
		int nj = j + direct[k][1];
		//out of boundary
		if (ni < 0 || nj < 0 || ni >= m.getSizeMaze() || nj >= m.rowSize()) continue;
		if (m.getMaze(ni, nj) == 9 || m.getMaze(ni, nj) == 10) continue;
		if (m.getMaze(ni, nj) == 0) count++;//visited
	}
	return count;
}

void MyMaze2dGenerator::rendomWalls(int a[], int size) {
	for (int i = 0; i < size; ++i)
	{
		//swap(a[i], a[rand() % size]);
		std::swap(visitOrderCell[i], visitOrderCell[rand() % size]);
	}
}


void MyMaze2dGenerator::generate_maze(int i, int j) {
	//set general Posituons thet is possible {up,dpwn,right,left}
	int direct[][2] = { {0,1}, {0,-1}, {-1,0}, {1,0} };

	//check out of the boundaries of the maze
	if (i < 0 || j < 0 || i >= m.getSizeMaze() || j >= m.rowSize()) return;

	//if visited, go back to the coming direction, return 
	if (m.getMaze(i, j) == 0) return;

	//If it is a cell initialized to start and end points then return to find another point to the wall
	if (m.getMaze(i, j) == 9 || m.getMaze(i, j) == 10) return;


	//some neightbors are visited in addition to the coming direction, return
	//this is to avoid circles in maze
	if (visitedAmount(i, j) > 1) return;
	if (m.getMaze(i, j) != 9 || m.getMaze(i, j) != 10)
		m.setMaze(i, j, 0);// set visited cube in the maze

	//shuffle the visitOrder
	rendomWalls(visitOrderCell,4); //mixed visitOrder Remdomly

	for (int k = 0; k < 4; ++k)
	{
		int ni = i + direct[visitOrderCell[k]][0];
		int nj = j + direct[visitOrderCell[k]][1];
		generate_maze(ni, nj);
	}
}

Maze2d& MyMaze2dGenerator::generate(int width, int height) {
	//Create a Maze in the Size the User requsted and set it as an object of the class!!
	Maze2d m11(width, height);
	m = m11;

	/*Terms for testing:
	**  Wall =1;
	**  Free = 0;
	**  PlayerPosition = 9;
	**  GoalPosition = 10;
	**  Visited = 0
	*/

	//First Set all maze into Walls Except from S and G
	int ii = 0;
	for (; ii < m.getSizeMaze(); ++ii)
	{
		for (int j = 0; j < m.rowSize(); ++j)
		{
			if (m.getMaze(ii, j) != 9 )
				if(m.getMaze(ii, j) != 10)
					m.setMaze(ii, j, 1); //Wall=1
		}
	}
	//DFS algoritem 
	//Realized in a separate function because it is recursive
	int i = 0;
	int j = 0;
	//send to function the first cube 
	//The function will go recursively over the whole matrix
	generate_maze(i, j);

	//Check if ther is aspace between S and G: if not so make one
	//TO DO //

	return m;
}

