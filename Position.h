#pragma once
#include<iostream> 
#include <stdlib.h>     /*	for rand metoed	 */
#include <vector>

			/********************************************/
			/************Part of task A.a****************/
			/***************Create Maze******************/
			/********************************************/

/*
** A class represents a movement of the user on the maze
*/

class Position
{
public:
	Position();
	Position* operator = (Position* p); //as copy constractor

	//Position* operator = (Maze2dSearchable* p);
	//Position(Maze2dSearchable& p);

	
	void choosePosition(int position);
	inline void StartPosition(int position) { _Startposition = position; }
	inline void CuurPosition(int position) { _Currposition = position; }
	inline void EndPosition(int position) { _Endposition = position; }
	inline int getXend() { return Xend; }
	inline int getYend() { return Yend; }


	inline int getPosition1() { return _position1; }
	inline int getPosition2() { return _position2; }
	inline void setX(int i) { _position1= i; }
	inline void setY(int j) { _position2 = j; }
	inline void setXend(int i) { Xend = i; }
	inline void setYend(int j) { Yend = j; }

	inline  const int getCuurPosition() { return _Currposition; }
	inline const int getStartPosition() { return _Startposition; }
	inline const int getEndPosition() { return _Endposition; }

	inline void setCurr(int curr) { _Currposition = curr; }
	inline void setStart(int start) { _Startposition = start; }
	inline void setEnd(int end) { _Endposition = end; }




	friend std::ostream& operator << (std::ostream& os, const Position& p) {
	//	os << " The Starting Point is in Cube :  " << p._Startposition << " " << std::endl;
		os << " The Current Point is in Cube : " << p._Currposition << " (x,y) = (" << p._position1 <<","<< p._position2 << ") " << std::endl;
		return os;
	}
//protected:   //if i do thet in the CommonSearcher it dosent work 
	int _position, _position1, _position2, _Startposition, _Currposition, _Endposition;
	int Xend, Yend;
};



