#pragma once
#include<iostream> 
#include <fstream>			 //ifstream
#include <map>			
#include<string>			/*  for string class  */ 




			/********************************************/
			/************Part of task B.b****************/
			/***************part of View*****************/
			/************Commendes from user*************/
			/********************************************/

class CLI
{
	CLI() {}
	CLI(std::ostream& out, std::istream& in); //in string from user 
	//friend std::ostream& operator << (std::ostream& os, const Maze2d& p) {}
	void start();

protected:
	std::map<int ,std::string> _data;			
};

