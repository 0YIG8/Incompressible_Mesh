#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <vector>
#include <string>
#include "Point.H"
#include "Face.H"

std::vector<Point> readPointsFile(std::string fileName);
std::vector<std::vector<int>> readFacesFile(std::string fileName);
std::vector<std::vector<int>> readCellsFile(std::string fileName);
void readBoundariesFile(std::string fileName, std::vector<Face>& listOfAllFaces);

#endif