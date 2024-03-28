#include "Point.H"

int Point::nPoints = 0; // Initialize no of pts to 0

/* 
    Constructor
*/

Point::Point(double x, double y, double z)
{
    this -> x = x;
    this -> y = y;
    this -> z = z;

    index = nPoints;
    nPoints++;
}



/* 
    Functions
*/

// Convert Point to Vector type
Vector Point::pointToVector() const
{
    Vector point(x, y, z);

    return point;
}



/* 
    Get Functions
*/

// Return point from index
Point Point::getPointFromIndex(int idx, const std::vector<Point> listOfAllPoints)
{
    Point p = listOfAllPoints[idx];
    return p;
}

int Point::getPointIndex() const
{
    return index;
}