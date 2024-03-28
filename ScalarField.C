#include "ScalarField.H"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>


/*
    Constructor
*/

ScalarField::ScalarField(std::vector<Cell>& listOfAllCells)
{
    size = listOfAllCells.size();  // Get number of cells in mesh
    hasBoundaryCondition = false;

    // Resize all matrices -> where, resize(nrow, ncol)
    coeffs.resize(size, size);
    source.resize(size, 1);
    variable.resize(size, 1);
    residual.resize(size, 1);

    // Initialize all variables to zero
    coeffs.setZero();
    source.setZero();
    variable.setOnes();

    residualNorm = 0.0;
}



/*
    Functions
*/

// Set initial condition for the field variable -- for same value in all cells
void ScalarField::setInitialCondition(double value)
{
    variable *= value;
}


// Reset coeff matrix and source matrix (eg. to be used at the end of a time step)
void ScalarField::resetMatrices()
{
    coeffs.setZero();
    source.setZero();
}


// Set a boundary condition for the field
void ScalarField::setBoundaryCondition(std::string patchName, std::string boundaryConditionType, double value)
{
    hasBoundaryCondition = true;

    nameOfBoundaryPatch.push_back(patchName);
    typeOfBoundaryCondition.push_back(boundaryConditionType);
    boundaryValue.push_back(value);
}


// Solve the system
void ScalarField::solveEquation()
{
    // variable = coeffs.lu().solve(source);
    // variable = coeffs.llt().solve(source);
    variable = coeffs.fullPivHouseholderQr().solve(source);
}


// Output data to a file for a single time step
void ScalarField::outputTimeStepDataToFile(std::string variableName, double timeStep)
{
    std::ostringstream oss;
    oss << variableName << "_" << timeStep << ".txt";
    std::ofstream output (oss.str());

    output << "cellIndex timeStep " << variableName << std::endl;

    for (int i = 0; i < size; i++)
    {
        output << i << " " << timeStep << " " << variable(i) << std::endl;
    }
}


// Calculate residual norm
void ScalarField::calculateResidualNorm()
{
    double sum = 0.0;

    for (int i = 0; i < size; i++)
    {
        sum += residual(i);
    }

    residualNorm = sum;
}



/* 
    Get Functions
*/

int ScalarField::getMatrixSize()
{
    return size;
}


bool ScalarField::getIfHasBoundaryCondition()
{
    return hasBoundaryCondition;
}


// For negative terms
void ScalarField::multiplyBy(double value, std::string matrix)
{
    if (matrix == "coeffs")
    {
        coeffs * value;
    }
    else if (matrix == "source")
    {
        source * value;
    }
    else if (matrix == "variable")
    {
        variable * value;
    }
}


// Get residual norm
double ScalarField::getResidualNorm()
{
    return residualNorm;
}