#include <iostream>
#include <set> 
#include <iterator> 
#include <stack>
#include <vector>
#include "CommonSearcher.h"

CommonSearcher::CommonSearcher() :m_evaluateNodes(0) {
	//row = { -1, 0, 0, 1 };
	row[0] = -1;
	row[1] = 0;
	row[2] = 0;
	row[3] = 1;

	//col = { 0, -1, 1, 0 };
	col[0] = 0;
	col[1] = -1;
	col[2] = 1;
	col[3] = 0;
}


bool operator<(const Position& p1, const Position& p2)
{
	return p1._Currposition < p2._Currposition;
}

void CommonSearcher::print_queue() {
	while (!m_openList.empty()) {
		std::cout << m_openList.top() << " ";
		m_openList.pop();
	}
	std::cout << '\n';
}

/////***------------------BFS----------------------------***///

bool  BFS::isValid(Maze2dSearchable& s , std::vector<std::vector<int>>  visited, int row, int col)
{

	return (row >= 0) && (row < s.getMaze().getW()) && (col >= 0) && (col < s.getMaze().getH())
		&& s.getMaze().getMaze(row, col) == 0 && !visited[row][col]; //visited =0 meens unvisted
}

bool  BFS::isValidEnd(Maze2dSearchable& s, std::vector<std::vector<int>>  visited, int row, int col)
{
	return (row >= 0) && (row < s.getMaze().getW()) && (col >= 0) && (col < s.getMaze().getH())
		&& s.getMaze().getMaze(row, col) == 10;
}

Solution BFS::search(Searchable& s) {
	Solution solve2;

	//Function for routing by game type
	//The implementation of the function is by type of game so it should not have been covered here at all.
	std::cout << "Function for routing by game type" << std::endl;

	return solve2;
}


void BFS::drowFormat(Maze2dSearchable& s, std::vector<std::vector<int>>  visited) {
	for (int i = 0; i <s.getMaze().getW() ; i++)
	{
		std::cout << "  ";
		for (int j = 0; j < s.getMaze().getH(); j++)
		{
			std::cout << " " << visited[i][j] << "  ";
		}
		std::cout << " ";
		std::cout << '\n';
	}
}


Solution BFS::search(Maze2dSearchable& s) {
	Solution solve; //Holds the array of solutions

	//get the (x,y) of end Position:
	int x = s.getpositions().getXend();
	int y = s.getpositions().getYend();
	GoalFound = false;


	//push it as a starting point of the path souletion
	m_openList.push(s.getpositions()); //look on position1 position2 as start 

	// initially all cells are unvisited
	std::vector<std::vector<int>>  _visited;
	for (int l = 0 ; l < s.getMaze().getW(); l++) {
		std::vector<int> row;
		for (int ll = 0; ll < s.getMaze().getH(); ll++) {
			row.push_back(0);
		}
		_visited.push_back(row);
	}

	while (!m_openList.empty())
	{
		// pop front node from queue and process it
		Position position = m_openList.top();
		m_openList.pop();

		// (i, j) represents current cell and dist stores its
		// minimum distance from the source
		int i = position.getPosition1(), j = position.getPosition2();
		//int i = node.x, j = node.y, dist = node.dist;

		// if destination is found 
		if (i == x && j == y)
		{
			int min_dist = m_evaluateNodes;
			break;
		}

		// check for all 4 possible movements from current cell
		for (int k = 0; k < 4; k++)
		{
			//when get to the goal 
			if (isValidEnd(s, _visited, i + row[k], j + col[k])) {
				prev.push_back(std::make_unique<Position>(s.getpositions()));
				s.setNewPosition(i + row[k], j + col[k]);
				prev.push_back(std::make_unique<Position>(s.getpositions()));
				GoalFound = true;
			}
		}

		if (GoalFound) break;

		//  enqueue each valid movement
		for (int k = 0; k < 4; k++)
		{

			// get all un visited nodes
			if (isValid(s, _visited , i + row[k], j + col[k]))
			{	
				// mark next cell as visited and enqueue it
				_visited[i + row[k]][j + col[k]] = 2;
				prev.push_back(std::make_unique<Position>(s.getpositions()));

				//s.getMaze().setMaze(i + row[k], j + col[k], 2);

				m_evaluateNodes = m_evaluateNodes + 1;
				//put the node to array of nabres
				s.setNewPosition(i + row[k], j + col[k]);
				m_openList.push(s.getpositions());
			}
		}
	}
	//if not found return empty
	if( !GoalFound ) 
		return solve;  //empty

	//if found set in a sulotion and retuen
	Position endP = s.getpositions();
	auto nextPtr = std::make_unique<Position>(endP); //end position


	auto nn = prev.size();
	int n = (int)nn;
		
	for (nextPtr; prev.size() ;) {

		if (n == 0) {
			solve.setSolution(std::move(nextPtr));
			return solve;
		}
		--n;
		std::unique_ptr<Position> nextPtr = std::move(prev.at(n));
		int x = nextPtr.get()->getPosition1();
		int y = nextPtr.get()->getPosition2();

		solve.setSolution(std::move(nextPtr));
		if(s.getMaze().getMaze(x,y) != 10 && s.getMaze().getMaze(x, y) != 9)
			s.getMaze().setMaze( x , y , 2 );
	}
}



///***------------------AStar----------------------------***///


// A structure to hold the neccesary parameters 


// A Utility Function to check whether given cell (row, col) 
// is a valid cell or not. 
bool AStar::isValid(int row, int col)
{
    // Returns true if row number and column number 
    // is in range 
    return (row >= 0) && (row < ROW) &&
        (col >= 0) && (col < COL);
}

// A Utility Function to check whether the given cell is 
// blocked or not 
bool AStar::isUnBlocked(Maze2dSearchable& s,int row, int col)
{
    // Returns true if the cell is not blocked else false 
    if (s[row][] == 1)
        return (true);
    else
        return (false);
}

// A Utility Function to check whether destination cell has 
// been reached or not 
bool AStar::isDestination(int row, int col, Cell::Pair dest)
{
    if (row == dest.first && col == dest.second)
        return (true);
    else
        return (false);
}

// A Utility Function to calculate the 'h' heuristics. 
double AStar::calculateHValue(int row, int col, Cell::Pair dest)
{
    // Return using the distance formula 
    return ((double)sqrt((row - dest.first) * (row - dest.first)
        + (col - dest.second) * (col - dest.second)));
}

// A Utility Function to trace the path from the source 
// to destination 
void AStar::tracePath(Cell cellDetails[][COL], Cell:: Pair dest)
{
    cout << "The path is: ";
    int row = dest.first;
    int col = dest.second;

    stack<Cell::Pair> Path;

    while (!(cellDetails[row][col].parent_i == row
        && cellDetails[row][col].parent_j == col))
    {
        Path.push(make_pair(row, col));
        int temp_row = cellDetails[row][col].parent_i;
        int temp_col = cellDetails[row][col].parent_j;
        row = temp_row;
        col = temp_col;
    }

    Path.push(make_pair(row, col));
    while (!Path.empty())
    {
        Cell::Pair p =  Path.top();
        Path.pop();
        cout << "->" << p.first << p.second;
    }

    return;
}

// A Function to find the shortest path between 
// a given source cell to a destination cell according 
// to A* Search Algorithm 
void AStar::aStarSearch(Cell cellDetails[][COL],Maze2dSearchable& s, Cell::Pair src, Cell::Pair dest)
{
    // If the source is out of range 
    if (isValid(src.first, src.second) == false)
    {
        
        cout << "source is invalid\n";
        return;
    }

    // If the destination is out of range 
    if (isValid(dest.first, dest.second) == false)
    {
        cout << "Destination is invalid\n";
        return;
    }

    // Either the source or the destination is blocked 
    if (isUnBlocked(s, src.first, src.second) == false ||
        isUnBlocked(s, dest.first, dest.second) == false)
    {
        cout << "Source or the destination is blocked\n";
        return;
    }

    // If the destination cell is the same as source cell 
    if (isDestination(src.first, src.second, dest) == true)
    {
       
        cout << "We are already at the destination\n";
        return;
    }

    // Create a closed list and initialise it to false which means 
    // that no cell has been included yet 
    // This closed list is implemented as a boolean 2D array 
    bool closedList[ROW][COL];
    memset(closedList, false, sizeof(closedList));

    // Declare a 2D array of structure to hold the details 
    //of that cell 
    //cell cellDetails[ROW][COL];

    int i, j;

    for (i = 0; i < *row; i++)
    {
        for (j = 0; j < *col; j++)
        {
            cellDetails[i][j].f = FLT_MAX;
            cellDetails[i][j].g = FLT_MAX;
            cellDetails[i][j].h = FLT_MAX;
            cellDetails[i][j].parent_i = -1;
            cellDetails[i][j].parent_j = -1;
        }
    }

    // Initialising the parameters of the starting node 
    i = src.first, j = src.second;
    cellDetails[i][j].f = 0.0;
    cellDetails[i][j].g = 0.0;
    cellDetails[i][j].h = 0.0;
    cellDetails[i][j].parent_i = i;
    cellDetails[i][j].parent_j = j;

    /*
     Create an open list having information as-
     <f, <i, j>>
     where f = g + h,
     and i, j are the row and column index of that cell
     Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
     This open list is implenented as a set of pair of pair.*/
    set<Cell::pPair> openList;

    // Put the starting cell on the open list and set its 
    // 'f' as 0 
    openList.insert(make_pair(0.0, make_pair(i, j)));

    // We set this boolean value as false as initially 
    // the destination is not reached. 
    bool foundDest = false;

    while (!openList.empty())
    {
        Cell::pPair p = *openList.begin();

        // Remove this vertex from the open list 
        openList.erase(openList.begin());

        // Add this vertex to the closed list 
        i = p.second.first;
        j = p.second.second;
        closedList[i][j] = true;

        /*
         Generating all the 8 successor of this cell

             N.W   N   N.E
               \   |   /
                \  |  /
             W----Cell----E
                  / | \
                /   |  \
             S.W    S   S.E

         Cell-->Popped Cell (i, j)
         N -->  North       (i-1, j)
         S -->  South       (i+1, j)
         E -->  East        (i, j+1)
         W -->  West           (i, j-1)
         N.E--> North-East  (i-1, j+1)
         N.W--> North-West  (i-1, j-1)
         S.E--> South-East  (i+1, j+1)
         S.W--> South-West  (i+1, j-1)*/

         // To store the 'g', 'h' and 'f' of the 8 successors 
        double gNew, hNew, fNew;

        //----------- 1st Successor (North) ------------ 
        
        // Only process this cell if this is a valid one 
        if (isValid(i - 1, j) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i - 1, j, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i - 1][j].parent_i = i;
                cellDetails[i - 1][j].parent_j = j;
                cout << "The destination cell is found\n";
                tracePath(cellDetails*, dest);
                foundDest = true;
                return;
            }
            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i - 1][j] == false &&
                isUnBlocked(s, i - 1, j) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i - 1, j, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i - 1][j].f == FLT_MAX ||
                    cellDetails[i - 1][j].f > fNew)
                {
                    openList.insert(make_pair(fNew,
                        make_pair(i - 1, j)));

                    // Update the details of this cell 
                    cellDetails[i - 1][j].f = fNew;
                    cellDetails[i - 1][j].g = gNew;
                    cellDetails[i - 1][j].h = hNew;
                    cellDetails[i - 1][j].parent_i = i;
                    cellDetails[i - 1][j].parent_j = j;
                }
            }
        }

        //----------- 2nd Successor (South) ------------ 

        // Only process this cell if this is a valid one 
        if (isValid(i + 1, j) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i + 1, j, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i + 1][j].parent_i = i;
                cellDetails[i + 1][j].parent_j = j;
                cout << "The destination cell is found\n";
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }
            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i + 1][j] == false &&
                isUnBlocked(s, i + 1, j) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i + 1, j, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i + 1][j].f == FLT_MAX ||
                    cellDetails[i + 1][j].f > fNew)
                {
                    openList.insert(make_pair(fNew, make_pair(i + 1, j)));
                    // Update the details of this cell 
                    cellDetails[i + 1][j].f = fNew;
                    cellDetails[i + 1][j].g = gNew;
                    cellDetails[i + 1][j].h = hNew;
                    cellDetails[i + 1][j].parent_i = i;
                    cellDetails[i + 1][j].parent_j = j;
                }
            }
        }

        //----------- 3rd Successor (East) ------------ 

        // Only process this cell if this is a valid one 
        if (isValid(i, j + 1) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i, j + 1, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i][j + 1].parent_i = i;
                cellDetails[i][j + 1].parent_j = j;
                cout << "The destination cell is found\n";
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i][j + 1] == false &&
                isUnBlocked(s, i, j + 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i, j + 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i][j + 1].f == FLT_MAX ||
                    cellDetails[i][j + 1].f > fNew)
                {
                    openList.insert(make_pair(fNew,
                        make_pair(i, j + 1)));

                    // Update the details of this cell 
                    cellDetails[i][j + 1].f = fNew;
                    cellDetails[i][j + 1].g = gNew;
                    cellDetails[i][j + 1].h = hNew;
                    cellDetails[i][j + 1].parent_i = i;
                    cellDetails[i][j + 1].parent_j = j;
                }
            }
        }

        //----------- 4th Successor (West) ------------ 

        // Only process this cell if this is a valid one 
        if (isValid(i, j - 1) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i, j - 1, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i][j - 1].parent_i = i;
                cellDetails[i][j - 1].parent_j = j;
               
                cout << "The destination cell found\n";
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i][j - 1] == false &&
                isUnBlocked(s, i, j - 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i, j - 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i][j - 1].f == FLT_MAX ||
                    cellDetails[i][j - 1].f > fNew)
                {
                    openList.insert(make_pair(fNew,
                        make_pair(i, j - 1)));

                    // Update the details of this cell 
                    cellDetails[i][j - 1].f = fNew;
                    cellDetails[i][j - 1].g = gNew;
                    cellDetails[i][j - 1].h = hNew;
                    cellDetails[i][j - 1].parent_i = i;
                    cellDetails[i][j - 1].parent_j = j;
                }
            }
        }

        //----------- 5th Successor (North-East) ------------ 

        // Only process this cell if this is a valid one 
        if (isValid(i - 1, j + 1) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i - 1, j + 1, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i - 1][j + 1].parent_i = i;
                cellDetails[i - 1][j + 1].parent_j = j;
               
                cout << "Destination cell found\n";
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i - 1][j + 1] == false &&
                isUnBlocked(s, i - 1, j + 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i - 1, j + 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i - 1][j + 1].f == FLT_MAX ||
                    cellDetails[i - 1][j + 1].f > fNew)
                {
                    openList.insert(make_pair(fNew,
                        make_pair(i - 1, j + 1)));

                    // Update the details of this cell 
                    cellDetails[i - 1][j + 1].f = fNew;
                    cellDetails[i - 1][j + 1].g = gNew;
                    cellDetails[i - 1][j + 1].h = hNew;
                    cellDetails[i - 1][j + 1].parent_i = i;
                    cellDetails[i - 1][j + 1].parent_j = j;
                }
            }
        }

        //----------- 6th Successor (North-West) ------------ 

        // Only process this cell if this is a valid one 
        if (isValid(i - 1, j - 1) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i - 1, j - 1, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i - 1][j - 1].parent_i = i;
                cellDetails[i - 1][j - 1].parent_j = j;
               
                cout << "Destination cell found\n";
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i - 1][j - 1] == false &&
                isUnBlocked(s, i - 1, j - 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i - 1, j - 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i - 1][j - 1].f == FLT_MAX ||
                    cellDetails[i - 1][j - 1].f > fNew)
                {
                    openList.insert(make_pair(fNew, make_pair(i - 1, j - 1)));
                    // Update the details of this cell 
                    cellDetails[i - 1][j - 1].f = fNew;
                    cellDetails[i - 1][j - 1].g = gNew;
                    cellDetails[i - 1][j - 1].h = hNew;
                    cellDetails[i - 1][j - 1].parent_i = i;
                    cellDetails[i - 1][j - 1].parent_j = j;
                }
            }
        }

        //----------- 7th Successor (South-East) ------------ 

        // Only process this cell if this is a valid one 
        if (isValid(i + 1, j + 1) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i + 1, j + 1, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i + 1][j + 1].parent_i = i;
                cellDetails[i + 1][j + 1].parent_j = j;
               
                cout << "Destination cell found\n";
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i + 1][j + 1] == false &&
                isUnBlocked(s, i + 1, j + 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i + 1, j + 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i + 1][j + 1].f == FLT_MAX ||
                    cellDetails[i + 1][j + 1].f > fNew)
                {
                    openList.insert(make_pair(fNew,
                        make_pair(i + 1, j + 1)));

                    // Update the details of this cell 
                    cellDetails[i + 1][j + 1].f = fNew;
                    cellDetails[i + 1][j + 1].g = gNew;
                    cellDetails[i + 1][j + 1].h = hNew;
                    cellDetails[i + 1][j + 1].parent_i = i;
                    cellDetails[i + 1][j + 1].parent_j = j;
                }
            }
        }

        //----------- 8th Successor (South-West) ------------ 

        // Only process this cell if this is a valid one 
        if (isValid(i + 1, j - 1) == true)
        {
            // If the destination cell is the same as the 
            // current successor 
            if (isDestination(i + 1, j - 1, dest) == true)
            {
                // Set the Parent of the destination cell 
                cellDetails[i + 1][j - 1].parent_i = i;
                cellDetails[i + 1][j - 1].parent_j = j;
                cout << "Destination cell found\n";
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed 
            // list or if it is blocked, then ignore it. 
            // Else do the following 
            else if (closedList[i + 1][j - 1] == false &&
                isUnBlocked(s, i + 1, j - 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i + 1, j - 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to 
                // the open list. Make the current square 
                // the parent of this square. Record the 
                // f, g, and h costs of the square cell 
                //                OR 
                // If it is on the open list already, check 
                // to see if this path to that square is better, 
                // using 'f' cost as the measure. 
                if (cellDetails[i + 1][j - 1].f == FLT_MAX ||
                    cellDetails[i + 1][j - 1].f > fNew)
                {
                    openList.insert(make_pair(fNew,
                        make_pair(i + 1, j - 1)));

                    // Update the details of this cell 
                    cellDetails[i + 1][j - 1].f = fNew;
                    cellDetails[i + 1][j - 1].g = gNew;
                    cellDetails[i + 1][j - 1].h = hNew;
                    cellDetails[i + 1][j - 1].parent_i = i;
                    cellDetails[i + 1][j - 1].parent_j = j;
                }
            }
        }
    }

    // When the destination cell is not found and the open 
    // list is empty, then we conclude that we failed to 
    // reach the destiantion cell. This may happen when the 
    // there is no way to destination cell (due to blockages) 
    if (foundDest == false)
        cout << "Failed to find destination\n";

    return;
}





