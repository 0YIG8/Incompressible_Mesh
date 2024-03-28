#include "Boundary.H"
#include "Face.H"
#include <string>


int Boundary::nBoundaryTypes = 0;

/* 
    Constructor
*/

Boundary::Boundary(std::vector<int> listOfBoundaryFaceIndices, std::string BCType, std::vector<Face>& listOfAllFaces)
{
    // Set as boundary faces
    for (size_t i = 0; i < listOfBoundaryFaceIndices.size(); i++)
    {
        for (size_t j = 0; j < listOfAllFaces.size(); j++)
        {
            if (listOfAllFaces[j].getFaceIndex() == listOfBoundaryFaceIndices[i])
            {
                listOfAllFaces[j].setAsBoundaryFace();
                listOfAllFaces[j].setBoundaryType(BCType);
            }
        }
    }

    nBoundaryTypes++;
}