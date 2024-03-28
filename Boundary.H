#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include <string>
#include "Face.H"


class Boundary
{
    static int nBoundaryTypes;

    private:
        std::vector<std::vector<int>> boundaryFaces;

    public:
        Boundary(std::vector<int> listOfBoundaryFaceIndices, std::string BCType, std::vector<Face>& listOfAllFaces);
};

#endif