#ifndef CELL_H
#define CELL_H

#include "Face.H"

class Cell
{
    static int nCells;  // Number of Cells

    private:
        std::vector<Face *> listOfFaces;

        double cellVolume;
        Vector cellCentre;
        int index;
        std::vector<int> listOfNeighbourCells;
        int numberOfCellNeighbours;

    public:
        Cell(std::vector<Face *> listOfFaces);

        Vector calcNaiveCellCentroid(std::vector<Face *> listOfFaces);
        void calcCellVolume();
        void calcCellCentre();
        void setFaceOwnerAndNeighbour(); // For faces of cells
        void setCellNeighbours(); 

        double getCellVolume() const;
        Vector getCellCentre() const;
        int getCellIndex() const;
        std::vector<int> getCellNeighbours() const;
        std::vector<Face *> getListOfFaces() const;
        int getNumberOfCellNeighbours() const;
};

#endif