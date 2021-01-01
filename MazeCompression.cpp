#include "MazeCompression.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

MazeCompression::MazeCompression(Maze2d& m) { //constarctor for empty maze
	my_maze = m.getData();  //Get data of maze and set it as object class
	w = m.getW();
	h = m.getH();
}

void MazeCompression::save(const char* FileName) {
	int xStart, yStart, xEnd, yEnd;
	//find the end and start position 
	for ( int i=0 ; i < w ; i++)
		for (int j = 0; j < h; j++) {
			if (my_maze[i][j] == 9) { //9 meens start position 
				xStart = i;
				yStart = j;
			}
			if (my_maze[i][j] == 10) { //109 meens goal position 
				xEnd = i;
				yEnd = j;
			}
		}


	ofstream myfile;
	myfile.open("example.txt", ios::trunc);  //open empty file for writing
	myfile << "SIZE: ( " << w << " , " << h << " )" << endl;
	myfile << "STAER: ( " << xStart << " , " << yStart << " )" << endl;
	myfile << "GOAL: ( " << xEnd << " , " << yEnd << " )" << endl;
	
	//Compression of the maze: 
	//markes : Wall=1 , Space =0 , Start = 9 , End =10 , PathSolution =2 (2 is not checked here)
	int count0 = 0, count1 = 0, findStart = 0, findEnd = 0; 
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if (my_maze[i][j] == 0)
				count0++;
			if (my_maze[i][j] == 1) 
				count0++;
			if (my_maze[i][j] == 9)
				count0++;
			if (my_maze[i][j] == 10)
				count0++;

		}
	}




	myfile.close();

}



void saveSolutionToStock();
Solution getExistSolution();

