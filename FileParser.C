#include "Vector.H"
#include "Point.H"
#include "Face.H"
#include "Boundary.H"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>


/* 
    Function that reads in a points file and returns a list of points 
*/
std::vector<Point> readPointsFile(std::string fileName)
{
    std::ifstream pointsFileIn;
    pointsFileIn.open(fileName);

    int nPts; // No. of points in the file
    std::string X, Y, Z; // To get and store x, y, z in string form
    double x, y, z; // To be passed to Point in vector
    std::string lineStr; // variable to save lines in string form
    std::vector<Point> listOfPoints;

    if (pointsFileIn.is_open() && !pointsFileIn.fail()) //if the file is open
    {
        std::getline(pointsFileIn, lineStr, '\n'); // Read no. of pts
        nPts = std::stoi(lineStr); // Store no of pts
        std::getline(pointsFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nPts; i++) //while the end of file is NOT reached
        {
            std::getline(pointsFileIn, lineStr, '(');
            std::getline(pointsFileIn, X, ' ');
            std::getline(pointsFileIn, Y, ' ');
            std::getline(pointsFileIn, Z, ')');

            x = std::stod(X);
            y = std::stod(Y);
            z = std::stod(Z);

            Point p(x, y, z);
            listOfPoints.push_back(p);
        }
    }
    else // If the file is not open
    {
        std::cout << "Unable to open " << fileName << "!" << std::endl;
    }

    return listOfPoints;
}



/* 
    Function that reads a list of indices of points in a face, for all faces, and returns it
*/
std::vector<std::vector<int>> readFacesFile(std::string fileName)
{
    std::ifstream facesFileIn;
    facesFileIn.open(fileName);

    int nFcs; // No. of faces in the file
    int nPts; // No. of points for each face
    std::string lineStr; // variable to save lines in string form
    std::vector<std::vector<int>> listOfFaces;

    if (facesFileIn.is_open()) //if the file is open
    {
        std::getline(facesFileIn, lineStr, '\n'); // Read no. of faces
        nFcs = std::stoi(lineStr); // Store no of faces
        std::getline(facesFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nFcs; i++)
        {   
            // Skip whitespace and read in number of points in the face
            std::getline(facesFileIn, lineStr, '(');
            std::istringstream ss(lineStr);
            ss >> std::ws;
            std::getline(ss, lineStr, '(');
            nPts = std::stoi(lineStr);

            std::vector<int> face; // Vector to store index of points for a single face
            int pointIdx; // Index of a point in a face

            // Store index of each point in the face std::vector
            for (int j = 0; j < nPts; j++)
            {
                if (j < nPts-1)
                {
                    std::getline(facesFileIn, lineStr, ' ');
                }
                else
                {
                    std::getline(facesFileIn, lineStr, '\n');
                }
                pointIdx = std::stoi(lineStr);
                face.push_back(pointIdx);
            }

            // Add the vector to the list of faces
            listOfFaces.push_back(face);
        }
    }
    else // If file is not open
    {
        std::cout << "Unable to open " << fileName << "!" << std::endl;
    }

    return listOfFaces;
}



/* 
    Function that reads a list of indices of faces in a cell, for all cells, and returns it 
*/
std::vector<std::vector<int>> readCellsFile(std::string fileName)
{
    std::ifstream cellsFileIn;
    cellsFileIn.open(fileName);

    int nCls; // No. of cells in the file
    int nFcs; // No. of faces belonging to the cell
    std::string lineStr; // variable to save lines in string form
    std::vector<std::vector<int>> listOfCells;

    if (cellsFileIn.is_open()) //if the file is open
    {
        std::getline(cellsFileIn, lineStr, '\n'); // Read no. of cells
        nCls = std::stoi(lineStr); // Store no of cells
        std::getline(cellsFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nCls; i++)
        {   
            // Skip whitespace and read in number of faces belonging to the cell
            std::getline(cellsFileIn, lineStr, '(');
            std::istringstream ss(lineStr);
            ss >> std::ws;
            std::getline(ss, lineStr, '(');
            nFcs = std::stoi(lineStr);

            std::vector<int> cell; // Vector to store index of faces for a single cell
            int faceIdx; // Index of a face in a cell

            // Store index of each face in the cell std::vector
            for (int j = 0; j < nFcs; j++)
            {
                if (j < nFcs-1)
                {
                    std::getline(cellsFileIn, lineStr, ' ');
                }
                else
                {
                    std::getline(cellsFileIn, lineStr, '\n');
                }
                faceIdx = std::stoi(lineStr);
                cell.push_back(faceIdx);
            }

            // Add the vector to the list of faces
            listOfCells.push_back(cell);
        }
    }
    else // If file is not open
    {
        std::cout << "Unable to open " << fileName << "!" << std::endl;
    }

    return listOfCells;
}



/* 
    Function that reads a list of boundary conditions and sets them
*/
void readBoundariesFile(std::string fileName, std::vector<Face>& listOfAllFaces)
{
    std::ifstream boundariesFileIn;
    boundariesFileIn.open(fileName);
    std::string lineStr; // variable to save lines in string form
    std::string bType;

    int nBCTypes = 0; // Number of differerent boundary conditions
    int nOfBFaces = 0; // Variable to store no. of boundary faces for each type 

    if (boundariesFileIn.is_open()) //if the file is open
    {
        std::getline(boundariesFileIn, lineStr, '\n'); // Read no. of BC types
        nBCTypes = std::stoi(lineStr); // Store no of BC types
        std::getline(boundariesFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nBCTypes; i++)
        {
            // Ignore whitespace and store boundary type
            std::getline(boundariesFileIn, bType, '\n');
            bType.erase(std::remove_if(bType.begin(), bType.end(), ::isspace), bType.end());

            std::getline(boundariesFileIn, lineStr, '\n');
            lineStr.erase(std::remove_if(lineStr.begin(), lineStr.end(), ::isspace), lineStr.end());
            nOfBFaces = std::stoi(lineStr);
            std::getline(boundariesFileIn, lineStr, '\n'); // Skip the line with the bracket

            // Get list of face indices of this boundary type
            std::vector<int> lf; // List of faces

            // Skip whitespace and read in indices of faces belonging to the cell
            std::getline(boundariesFileIn, lineStr, '\n');
            auto whitespace = std::find_if(lineStr.begin(), lineStr.end(), [](unsigned char ch) {return !std::isspace(ch);});
            lineStr.erase(lineStr.begin(), whitespace);
            std::istringstream ss(lineStr);
                
            for (int j = 0; j < nOfBFaces; j++)
            {
                if (j < nOfBFaces-1)
                {
                    std::getline(ss, lineStr, ' ');
                }
                else
                {
                    std::getline(ss, lineStr, '\n');
                }

                lf.push_back(std::stoi(lineStr));
            }

            std::getline(boundariesFileIn, lineStr, '\n'); // Skip the line with the bracket

            // Set up boundaries
            Boundary b(lf, bType, listOfAllFaces);
        }
    }
    else // If file is not open
    {
        std::cout << "Unable to open " << fileName << "!" << std::endl;
    }
}