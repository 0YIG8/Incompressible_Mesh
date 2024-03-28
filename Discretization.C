#include <iostream>
#include "Vector.H"
#include "Discretization.H"
#include "Face.H"
#include "Cell.H"
#include <Eigen/Dense>


/*
        DISCRETIZATION ESSENTIALS
*/
// Calc and return dist d between P and N (d = delta)
double delta(Face* f, std::vector<Cell>& listOfAllCells)
{
    double delta = 0.0;
    Vector ownerCentre = listOfAllCells[f->getFaceOwner()].getCellCentre();
    Vector faceCentre = f->getFaceCentre();

    if (f->getIfBoundaryFace() == true)
    {
        delta = (faceCentre - ownerCentre).norm();
    }
    else 
    {
        Vector neighbourCentre = listOfAllCells[f->getFaceNeighbour()].getCellCentre();
        delta = (faceCentre - ownerCentre).norm() + (faceCentre - neighbourCentre).norm();
    }

    return delta;
}


// Calc and return interpolation factor f_x
double interpolationFactor(Face* f, std::vector<Cell>& listOfAllCells)
{
    double interpolationFactor = 0.0;
    Vector faceCentre = f->getFaceCentre();

    if (f->getIfBoundaryFace() == true)
    {
        interpolationFactor = 1.0;
    }
    else
    {
        Vector neighbourCentre = listOfAllCells[f->getFaceNeighbour()].getCellCentre();
        interpolationFactor = (faceCentre - neighbourCentre).norm() / delta(f, listOfAllCells);
    }

    return interpolationFactor;
}


// Calculate Time Derivative Term  (For Scalar Field)
void evaluateTimeDerivative(std::vector<Cell>& listOfAllCells, ScalarField& fieldVariable, double delta_t)
{
    for (size_t i = 0; i < listOfAllCells.size(); i++)
    {
        Cell c = listOfAllCells[i];

        fieldVariable.coeffs(i, i) += c.getCellVolume() / delta_t;  // Diagonal contribution
        fieldVariable.source(i) += (c.getCellVolume() * fieldVariable.variable(i)) / delta_t;  // Source contribution
    }
}


// Calculate Time Derivative Term  (For Vector Field) --- NOT TESTED!!!
void evaluateTimeDerivative(std::vector<Cell>& listOfAllCells, VectorField& fieldVariable, double delta_t)
{
    for (size_t i = 0; i < listOfAllCells.size(); i++)
    {
        Cell c = listOfAllCells[i];

        fieldVariable.coeffs(i, i) += c.getCellVolume() / delta_t;  // Diagonal contribution

        fieldVariable.sourceX(i) += (c.getCellVolume() * fieldVariable.variableX(i)) / delta_t;  // Source contribution in X direction
        fieldVariable.sourceY(i) += (c.getCellVolume() * fieldVariable.variableY(i)) / delta_t;  // Source contribution in Y direction
        fieldVariable.sourceZ(i) += (c.getCellVolume() * fieldVariable.variableZ(i)) / delta_t;  // Source contribution in Z direction
    }
}


// Calculate Diffusion Term (For Scalar Field --- with constant gamma)
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, double gamma)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.source(idxP) += gamma * fieldVariable.boundaryValue[k] * f->getFaceArea() / f->getFaceDelta();
                        fieldVariable.coeffs(idxP, idxP) += -(-gamma * f->getFaceArea() / f->getFaceDelta());
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann")
                    {
                        fieldVariable.source(idxP) += gamma * f->getFaceArea() * fieldVariable.boundaryValue[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            double coeffP = gamma * f->getFaceArea() / f->getFaceDelta();
            double coeffN = -coeffP;

            fieldVariable.coeffs(idxP, idxP) += coeffP;
            fieldVariable.coeffs(idxP, idxN) += coeffN;
            fieldVariable.coeffs(idxN, idxP) += coeffN;
            fieldVariable.coeffs(idxN, idxN) += coeffP;
        }
    }
}


// Calculate Diffusion Term (For Scalar Field --- interpolated gamma from P and N cells)  --- for negative diffusion term
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, std::vector<double>& gamma)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();
        double fx = f->getFaceInterpolationFactor();

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            double gamma_interpolated = fx * gamma[idxP];

            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.source(idxP) += gamma_interpolated * fieldVariable.boundaryValue[k] * f->getFaceArea() / f->getFaceDelta();
                        fieldVariable.coeffs(idxP, idxP) += -(-gamma_interpolated * f->getFaceArea() / f->getFaceDelta());
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann")
                    {
                        fieldVariable.source(idxP) += gamma_interpolated * f->getFaceArea() * fieldVariable.boundaryValue[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            double gamma_interpolated = (fx * gamma[idxP]) + ((1.0 - fx) * gamma[idxN]);

            double coeffP = gamma_interpolated * f->getFaceArea() / f->getFaceDelta();
            double coeffN = -coeffP;

            fieldVariable.coeffs(idxP, idxP) += coeffP;
            fieldVariable.coeffs(idxP, idxN) += coeffN;
            fieldVariable.coeffs(idxN, idxP) += coeffN;
            fieldVariable.coeffs(idxN, idxN) += coeffP;
        }
    }
}



// Calculate Diffusion Term (For Scalar Field --- interpolated gamma from P and N cells)  --- for positive diffusion term
void evaluatePositiveDiffusionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, std::vector<double>& gamma)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();
        double fx = f->getFaceInterpolationFactor();

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            double gamma_interpolated = fx * gamma[idxP];

            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.source(idxP) += -gamma_interpolated * fieldVariable.boundaryValue[k] * f->getFaceArea() / f->getFaceDelta();
                        fieldVariable.coeffs(idxP, idxP) += (-gamma_interpolated * f->getFaceArea() / f->getFaceDelta());
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann")
                    {
                        fieldVariable.source(idxP) += -gamma_interpolated * f->getFaceArea() * fieldVariable.boundaryValue[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            double gamma_interpolated = (fx * gamma[idxP]) + ((1.0 - fx) * gamma[idxN]);

            double coeffN = gamma_interpolated * f->getFaceArea() / f->getFaceDelta();
            double coeffP = -coeffN;

            fieldVariable.coeffs(idxP, idxP) += coeffP;
            fieldVariable.coeffs(idxP, idxN) += coeffN;
            fieldVariable.coeffs(idxN, idxP) += coeffN;
            fieldVariable.coeffs(idxN, idxN) += coeffP;
        }
    }
}



// Calculate Diffusion Term (For Vector Field --- with constant gamma)
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, double gamma)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.sourceX(idxP) += gamma * fieldVariable.boundaryValueX[k] * f->getFaceArea() / f->getFaceDelta();
                        fieldVariable.sourceY(idxP) += gamma * fieldVariable.boundaryValueY[k] * f->getFaceArea() / f->getFaceDelta();
                        fieldVariable.sourceZ(idxP) += gamma * fieldVariable.boundaryValueZ[k] * f->getFaceArea() / f->getFaceDelta();

                        fieldVariable.coeffs(idxP, idxP) += (gamma * f->getFaceArea() / f->getFaceDelta());
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann")
                    {
                        fieldVariable.sourceX(idxP) += gamma * f->getFaceArea() * fieldVariable.boundaryValueX[k];
                        fieldVariable.sourceY(idxP) += gamma * f->getFaceArea() * fieldVariable.boundaryValueY[k];
                        fieldVariable.sourceZ(idxP) += gamma * f->getFaceArea() * fieldVariable.boundaryValueZ[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            double coeffP = gamma * f->getFaceArea() / f->getFaceDelta();
            double coeffN = -coeffP;

            fieldVariable.coeffs(idxP, idxP) += coeffP;
            fieldVariable.coeffs(idxP, idxN) += coeffN;
            fieldVariable.coeffs(idxN, idxP) += coeffN;
            fieldVariable.coeffs(idxN, idxN) += coeffP;
        }
    }
}


// Calculate Diffusion Term (For Vector Field --- with different gamma for every cell)
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, std::vector<double>& gamma)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.sourceX(idxP) += gamma[idxP] * fieldVariable.boundaryValueX[k] * f->getFaceArea() / f->getFaceDelta();
                        fieldVariable.sourceY(idxP) += gamma[idxP] * fieldVariable.boundaryValueY[k] * f->getFaceArea() / f->getFaceDelta();
                        fieldVariable.sourceZ(idxP) += gamma[idxP] * fieldVariable.boundaryValueZ[k] * f->getFaceArea() / f->getFaceDelta();

                        fieldVariable.coeffs(idxP, idxP) += -(-gamma[idxP] * f->getFaceArea() / f->getFaceDelta());
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann")
                    {
                        fieldVariable.sourceX(idxP) += gamma[idxP] * f->getFaceArea() * fieldVariable.boundaryValueX[k];
                        fieldVariable.sourceY(idxP) += gamma[idxP] * f->getFaceArea() * fieldVariable.boundaryValueY[k];
                        fieldVariable.sourceZ(idxP) += gamma[idxP] * f->getFaceArea() * fieldVariable.boundaryValueZ[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            double fx = f->getFaceInterpolationFactor();
            double gamma_interpolated = (fx * gamma[idxP]) + (1.0 - fx) * gamma[idxN];

            double coeffP = gamma_interpolated * f->getFaceArea() / f->getFaceDelta();
            double coeffN = -coeffP;

            fieldVariable.coeffs(idxP, idxP) += coeffP; 
            fieldVariable.coeffs(idxP, idxN) += coeffN; 
            fieldVariable.coeffs(idxN, idxP) += coeffN; 
            fieldVariable.coeffs(idxN, idxN) += coeffP; 
        }
    }
}


// Calculate Convection Term (For Scalar Field)
void evaluateConvectionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, Vector& constVelocity, std::string scheme)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();

        double F = getFlux(f, constVelocity);
        f->setFaceFlux(F);

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.source(idxP) += -F * fieldVariable.boundaryValue[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann")
                    {
                        fieldVariable.source(idxP) += -F * f->getFaceDelta() * fieldVariable.boundaryValue[k];
                        fieldVariable.coeffs(idxP, idxP) += F;
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            if (scheme == "CD" || scheme == "Central Differencing")
            {
                double coeffP = f->getFaceInterpolationFactor() * F;
                double coeffN = (1.0 - f->getFaceInterpolationFactor()) * F;

                fieldVariable.coeffs(idxP, idxP) += coeffP;
                fieldVariable.coeffs(idxP, idxN) += coeffN;
                fieldVariable.coeffs(idxN, idxP) += -coeffP;
                fieldVariable.coeffs(idxN, idxN) += -coeffN;
            }
            else if (scheme == "UD" || scheme == "Upwind Differencing" || scheme == "Upwind")
            {
                double coeffP = pos(F) * f->getFaceFlux();
                double coeffN = (1.0 - pos(F)) * f->getFaceFlux();

                fieldVariable.coeffs(idxP, idxP) += coeffP;
                fieldVariable.coeffs(idxP, idxN) += coeffN;
                fieldVariable.coeffs(idxN, idxP) += -coeffP;
                fieldVariable.coeffs(idxN, idxN) += -coeffN;
            }
            else // If none of these
            {
                std::cout << "Invalid convection scheme!" << std::endl;
                exit(1);
            }
        }
    }
}


// Calculate Convection Term (For Vector Field + constant velocity)
void evaluateConvectionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, Vector& constVelocity, std::string scheme)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();

        double F = getFlux(f, constVelocity);
        f->setFaceFlux(F);

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.sourceX(idxP) += -F * fieldVariable.boundaryValueX[k];
                        fieldVariable.sourceY(idxP) += -F * fieldVariable.boundaryValueY[k];
                        fieldVariable.sourceZ(idxP) += -F * fieldVariable.boundaryValueZ[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann")
                    {
                        fieldVariable.sourceX(idxP) += -F * f->getFaceDelta() * fieldVariable.boundaryValueX[k];
                        fieldVariable.sourceY(idxP) += -F * f->getFaceDelta() * fieldVariable.boundaryValueY[k];
                        fieldVariable.sourceZ(idxP) += -F * f->getFaceDelta() * fieldVariable.boundaryValueZ[k];

                        fieldVariable.coeffs(idxP, idxP) += F;
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            if (scheme == "CD" || scheme == "Central Differencing")
            {
                double coeffP = f->getFaceInterpolationFactor() * F;
                double coeffN = (1.0 - f->getFaceInterpolationFactor()) * F;

                fieldVariable.coeffs(idxP, idxP) += coeffP;
                fieldVariable.coeffs(idxP, idxN) += coeffN;
                fieldVariable.coeffs(idxN, idxP) += -coeffP;
                fieldVariable.coeffs(idxN, idxN) += -coeffN;
            }
            else if (scheme == "UD" || scheme == "Upwind Differencing" || scheme == "Upwind")
            {
                double coeffP = pos(F) * F;
                double coeffN = (1.0 - pos(F)) * F;

                fieldVariable.coeffs(idxP, idxP) += coeffP;
                fieldVariable.coeffs(idxP, idxN) += coeffN;
                fieldVariable.coeffs(idxN, idxP) += -coeffP;
                fieldVariable.coeffs(idxN, idxN) += -coeffN;
            }
            else // If none of these
            {
                std::cout << "Invalid convection scheme!" << std::endl;
                exit(1);
            }
        }
    }
}


// Calculate Convection Term (For Vector Field + velocity field)
void evaluateConvectionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, VectorField& velocityField, std::string scheme)
{
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();

        double F = getFlux(f, velocityField, idxP, idxN);
        double Flux = f->getFaceFlux();

        if (f->getIfBoundaryFace() == true) // For boundary faces
        {
            for (size_t k = 0; k < fieldVariable.nameOfBoundaryPatch.size(); k++)
            {
                if (f->getBoundaryType() == fieldVariable.nameOfBoundaryPatch[k])
                {
                    if (fieldVariable.typeOfBoundaryCondition[k] == "Dirichlet")
                    {
                        fieldVariable.sourceX(idxP) += -Flux * fieldVariable.boundaryValueX[k];  
                        fieldVariable.sourceY(idxP) += -Flux * fieldVariable.boundaryValueY[k];
                        fieldVariable.sourceZ(idxP) += -Flux * fieldVariable.boundaryValueZ[k];
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Neumann") 
                    {
                        fieldVariable.sourceX(idxP) += -Flux * f->getFaceDelta() * fieldVariable.boundaryValueX[k];
                        fieldVariable.sourceY(idxP) += -Flux * f->getFaceDelta() * fieldVariable.boundaryValueY[k];
                        fieldVariable.sourceZ(idxP) += -Flux * f->getFaceDelta() * fieldVariable.boundaryValueZ[k];

                        fieldVariable.coeffs(idxP, idxP) += Flux;
                    }
                    else if (fieldVariable.typeOfBoundaryCondition[k] == "Mixed")
                    {
                        // Do something
                    }
                    else if (f->getBoundaryType() == "empty")
                    {
                        break;
                    }
                    else // Not any recognised boundary condition
                    {
                        std::cout << "Invalid boundary condition!" << std::endl;
                    }
                }
            }
        }
        else // If not a boundary face
        {
            if (scheme == "CD" || scheme == "Central Differencing")
            {
                double coeffP = f->getFaceInterpolationFactor() * Flux;
                double coeffN = (1.0 - f->getFaceInterpolationFactor()) * Flux;

                fieldVariable.coeffs(idxP, idxP) += coeffP;  
                fieldVariable.coeffs(idxP, idxN) += coeffN; 
                fieldVariable.coeffs(idxN, idxP) += -coeffP;
                fieldVariable.coeffs(idxN, idxN) += -coeffN;
            }
            else if (scheme == "UD" || scheme == "Upwind Differencing" || scheme == "Upwind")
            {
                double coeffP = pos(f->getFaceFlux()) * Flux;
                double coeffN = (1.0 - pos(f->getFaceFlux())) * Flux;

                fieldVariable.coeffs(idxP, idxP) += coeffP;
                fieldVariable.coeffs(idxP, idxN) += coeffN;
                fieldVariable.coeffs(idxN, idxP) += -coeffP;
                fieldVariable.coeffs(idxN, idxN) += -coeffN;
            }
            else // If none of these
            {
                std::cout << "Invalid convection scheme!" << std::endl;
                exit(1);
            }
        }
    }
}


/*
        FLUX CALCULATION
*/
double getFlux(Face* f, Vector constVelocity) // For const velocity
{
    double flux = (f->getFaceNormal() * f->getFaceArea()).dotMult(constVelocity);  // F = sf.dot_mult(velocity)

    return flux;
}

// Specific for dirichlet condition right now
double getFlux(Face* f, VectorField& velocityField, int& idxP, int& idxN) // For vector field velocity (interpolated u)
{
    double fx = f->getFaceInterpolationFactor();
    double vel_X = 0.0, vel_Y = 0.0, vel_Z = 0.0;

    // Interpolate velocity from owner and neighbour of face
    if (f->getIfBoundaryFace() == true)  // If the face is a boundary
    {
        for (size_t k = 0; k < velocityField.nameOfBoundaryPatch.size(); k++)
        {
            if (f->getBoundaryType() == velocityField.nameOfBoundaryPatch[k])
            {
                if (velocityField.typeOfBoundaryCondition[k] == "Dirichlet")
                {
                    vel_X = velocityField.boundaryValueX[k];  
                    vel_Y = velocityField.boundaryValueY[k];
                    vel_Z = velocityField.boundaryValueZ[k];
                }
                else if (velocityField.typeOfBoundaryCondition[k] == "Neumann")   
                {
                    // Do something
                }
                else if (velocityField.typeOfBoundaryCondition[k] == "Mixed")
                {
                    // Do something
                }
                else if (f->getBoundaryType() == "empty")
                {
                    break;
                }
                else // Not any recognised boundary condition
                {
                    std::cout << "Invalid boundary condition!" << std::endl;
                }
            }
        }
    }
    else  // If not a boundary face
    {
        vel_X = (fx * velocityField.variableX(idxP)) + ((1.0 - fx) * velocityField.variableX(idxN));
        vel_Y = (fx * velocityField.variableY(idxP)) + ((1.0 - fx) * velocityField.variableY(idxN));
        vel_Z = (fx * velocityField.variableZ(idxP)) + ((1.0 - fx) * velocityField.variableZ(idxN));
    }

    Vector vel(vel_X, vel_Y, vel_Z);

    // Calculate flux
    double flux = (f->getFaceNormal() * f->getFaceArea()).dotMult(vel);

    return flux;
}


double pos(double flux)  // For upwind differencing 
{
    if (flux >= 0.0)
    {
        return 1.0;
    }
    else // if (flux < 0.0)
    {
        return 0.0;
    }
}



/*
        GRADIENT CALCULATION
*/
// Calc and return Gradient for a cell using Gauss Theorem 
Vector getGradientGauss(ScalarField& fieldVariable, int cellIdx, std::vector<Cell>& listOfAllCells)
{
    Vector gradient(0.0, 0.0, 0.0);

    Cell c = listOfAllCells[cellIdx];
    std::vector<Face *> lf = c.getListOfFaces();  // Get faces for this cell

    for (int i = 0; i < lf.size(); i++)
    {
        Face* f = lf[i];  // Pick a face from the list
        int idxP = f->getFaceOwner();
        int idxN = f->getFaceNeighbour();
        double fx = f->getFaceInterpolationFactor();

        if (f->getIfBoundaryFace() == true) // if boundary face
        {
            gradient = gradient + (f->getFaceNormal() * f->getFaceArea() * (fx * fieldVariable.variable(idxP)));
        }
        else // if not boundary face
        {
            if(cellIdx == idxP)
            {
                gradient = gradient + (f->getFaceNormal() * f->getFaceArea() * ((fx * fieldVariable.variable(idxP)) + ((1.0 - fx) * fieldVariable.variable(idxN))));
            }
            else
            {
                gradient = gradient - (f->getFaceNormal() * f->getFaceArea() * ((fx * fieldVariable.variable(idxP)) + ((1.0 - fx) * fieldVariable.variable(idxN))));
            }
            
        }
    }

    return gradient;
}



/* 
        UNDER RELAXATION
*/
void underRelaxImplicit(VectorField& fieldVariable, double alpha)
{
    int size = fieldVariable.getMatrixSize();

    for (int i = 0; i < size; i++)
    {
        double a_P = fieldVariable.coeffs(i, i);

        fieldVariable.sourceX(i) += ((1.0 - alpha) / alpha) * a_P * fieldVariable.variableX(i);
        fieldVariable.sourceY(i) += ((1.0 - alpha) / alpha) * a_P * fieldVariable.variableY(i);
        fieldVariable.sourceZ(i) += ((1.0 - alpha) / alpha) * a_P * fieldVariable.variableZ(i);

        // fieldVariable.coeffs(i, i) /= alpha;
        fieldVariable.coeffs(i, i) = fieldVariable.coeffs(i, i) / alpha;
    }
}


void underRelaxExplicit(ScalarField& fieldVariable, Eigen::VectorXd& p_old, double alpha)
{
    int size = fieldVariable.getMatrixSize();

    for (int i = 0; i < size; i++)
    {
        fieldVariable.variable(i) = p_old(i) + (alpha * (fieldVariable.variable(i) - p_old(i)));   // p** = p* + alpha(p - p*)
    }
}


/*
    CALCULATE RESIDUALS
*/
void calculateResidual(VectorField fieldVariable)
{
    for (int i = 0; i < fieldVariable.getMatrixSize(); i++)
    {
        fieldVariable.residualX = fieldVariable.sourceX - (fieldVariable.coeffs * fieldVariable.variableX);
        fieldVariable.residualY = fieldVariable.sourceY - (fieldVariable.coeffs * fieldVariable.variableY);
        fieldVariable.residualZ = fieldVariable.sourceZ - (fieldVariable.coeffs * fieldVariable.variableZ);
    }
}

void calculateResidual(ScalarField fieldVariable)
{
    for (int i = 0; i < fieldVariable.getMatrixSize(); i++)
    {
        fieldVariable.residual = fieldVariable.source - (fieldVariable.coeffs * fieldVariable.variable);
    }
}