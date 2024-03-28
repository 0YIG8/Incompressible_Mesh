#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <vector>
#include "Cell.H"
#include <Eigen/Dense>

class VectorField
{
    private:
        int size;
        bool hasBoundaryCondition;

    public:
        VectorField(std::vector<Cell>& listOfAllCells);

        Eigen::MatrixXd coeffs; // Matrix of dimensions dynamic x dynamic of type double
        Eigen::VectorXd variableX; // Matrix of dimensions dynamic x 1 of type double
        Eigen::VectorXd variableY;
        Eigen::VectorXd variableZ;
        Eigen::VectorXd sourceX;
        Eigen::VectorXd sourceY;
        Eigen::VectorXd sourceZ;

        Eigen::VectorXd residualX;
        Eigen::VectorXd residualY;
        Eigen::VectorXd residualZ;
        std::vector<double> residualNorm;

        std::vector<std::string> nameOfBoundaryPatch;
        std::vector<std::string> typeOfBoundaryCondition;
        std::vector<double> boundaryValueX;
        std::vector<double> boundaryValueY;
        std::vector<double> boundaryValueZ;

        void setInitialCondition(double valueX, double valueY, double valueZ);
        void resetMatrices();
        void setBoundaryCondition(std::string patchName, std::string boundaryConditionType, double valueX, double valueY, double valueZ);
        void solveEquation();

        void calculateResidualNorm();

        int getMatrixSize();
        bool getIfHasBoundaryCondition();
        std::vector<double> getResidualNorm();
};

#endif