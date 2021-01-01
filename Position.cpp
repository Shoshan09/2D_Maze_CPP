#include "Position.h"


Position::Position() : _position(0), _position1(0), _position2(0), _Startposition(0), _Currposition(0), _Endposition(0){}

Position* Position::operator = (Position* p){
	_Currposition = p->_Currposition;
	_Startposition = p->_Startposition;
	_Endposition = p->_Endposition;
	_position1 = p->_position1;
	_position2 = p->_position2;
	
	return this;
}

//Position::Position(Maze2dSearchable& p) {
//	int x;
//	//p.getPositions().setCurr(_Currposition);
//	//p.getPositions().setStart(_Startposition);
//	//p.getPositions().setEnd(_Endposition);
//	//p.getPositions().setX(_position1);
//	//p.getPositions().setY(_position2);
//}


void Position::choosePosition(int position)
{
	_position = position;
}


