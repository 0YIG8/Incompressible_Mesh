#include "VectorField.H"
#include <Eigen/Dense>


/*
    Constructor
*/

VectorField::VectorField(std::vector<Cell>& listOfAllCells)
{
    size = listOfAllCells.size();  // Get number of cells in mesh

    // Resize all matrices -> where, resize(nrow, ncol)
    coeffs.resize(size, size);
    sourceX.resize(size, 1);
    sourceY.resize(size, 1);
    sourceZ.resize(size, 1);
    variableX.resize(size, 1);
    variableY.resize(size, 1);
    variableZ.resize(size, 1);
    residualX.resize(size, 1);
    residualY.resize(size, 1);
    residualZ.resize(size, 1);
    residualNorm.resize(3);

    // Initialize all variables to zero
    coeffs.setZero();
    sourceX.setZero();
    sourceY.setZero();
    sourceZ.setZero();
    variableX.setOnes();
    variableY.setOnes();
    variableZ.setOnes();
}



/*
    Functions
*/

// Set initial condition for the field variable -- for same value in all cells
void VectorField::setInitialCondition(double valueX, double valueY, double valueZ)
{
    variableX *= valueX;
    variableY *= valueY;
    variableZ *= valueZ;
}


// Reset coeff matrix and source matrix (eg. to be used at the end of a time step)
void VectorField::resetMatrices()
{
    coeffs.setZero();

    sourceX.setZero();
    sourceY.setZero();
    sourceZ.setZero();
}


// Set a boundary condition for the field
void VectorField::setBoundaryCondition(std::string patchName, std::string boundaryConditionType, double valueX, double valueY, double valueZ)
{
    hasBoundaryCondition = true;

    nameOfBoundaryPatch.push_back(patchName);
    typeOfBoundaryCondition.push_back(boundaryConditionType);
    boundaryValueX.push_back(valueX);
    boundaryValueY.push_back(valueY);
    boundaryValueZ.push_back(valueZ);
}


// Solve the system
void VectorField::solveEquation()
{
    variableX = coeffs.lu().solve(sourceX);
    variableY = coeffs.lu().solve(sourceY);
    variableZ = coeffs.lu().solve(sourceZ);
}


// Calculate residual norm
void VectorField::calculateResidualNorm()
{
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;

    for (int i = 0; i < size; i++)
    {
        sumX += residualX(i);
        sumY += residualY(i);
        sumZ += residualZ(i);
    }

    residualNorm[0] = sumX;
    residualNorm[1] = sumY;
    residualNorm[2] = sumZ;
}


/* 
    Get Functions
*/

int VectorField::getMatrixSize()
{
    return size;
}


bool VectorField::getIfHasBoundaryCondition()
{
    return hasBoundaryCondition;
}

std::vector<double> VectorField::getResidualNorm()
{
    return residualNorm;
}