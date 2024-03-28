/* 
    P-U COUPLING FOR THE INCOMPRESSIBLE NAVIER STOKES EQUATIONS USING THE SIMPLE ALGORITHM
    By BCN - 0YIG8
*/

#include "Vector.H"
#include "Mesh.H"
#include "Discretization.H"
#include "ScalarField.H"
#include "OutputData.H"
#include <ctime>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>


int main(void)
{
    /*
        RECORD SIMULATION TIME
    */
    clock_t tStart, tEnd;
    tStart = clock();


    /*
        SPECIFY INPUT FILES FOR MESH
    */
    // std::string PointsFile = "points_5x5.txt";
    // std::string FacesFile = "faces_5x5.txt";
    // std::string CellsFile = "cells_5x5.txt";
    // std::string BoundariesFile = "boundary_5x5.txt";  
    // int nCellsx = 5, nCellsy = 5;

    std::string PointsFile = "40x40/40x40_points.txt";
    std::string FacesFile = "40x40/40x40_faces.txt";
    std::string CellsFile = "40x40/40x40_cells.txt";
    std::string BoundariesFile = "40x40/40x40_boundaries.txt";
    int nCellsx = 40, nCellsy = 40;

    // std::string PointsFile = "60x60/60x60_points.txt";
    // std::string FacesFile = "60x60/60x60_faces.txt";
    // std::string CellsFile = "60x60/60x60_cells.txt";
    // std::string BoundariesFile = "60x60/60x60_boundaries.txt";
    // int nCellsx = 60, nCellsy = 60;

    // std::string PointsFile = "80x80/80x80_points.txt";
    // std::string FacesFile = "80x80/80x80_faces.txt";
    // std::string CellsFile = "80x80/80x80_cells.txt";
    // std::string BoundariesFile = "80x80/80x80_boundaries.txt";
    // int nCellsx = 80, nCellsy = 80;

    // std::string PointsFile = "120x120/120x120_points.txt";
    // std::string FacesFile = "120x120/120x120_faces.txt";
    // std::string CellsFile = "120x120/120x120_cells.txt";
    // std::string BoundariesFile = "120x120/120x120_boundaries.txt";
    // int nCellsx = 120, nCellsy = 120;

    // std::string PointsFile = "160x160/160x160_points.txt";
    // std::string FacesFile = "160x160/160x160_faces.txt";
    // std::string CellsFile = "160x160/160x160_cells.txt";
    // std::string BoundariesFile = "160x160/160x160_boundaries.txt";
    // int nCellsx = 160, nCellsy = 160;



    /*
        READ IN FILES TO GENERATE THE MESH
    */
    Mesh cavity_mesh(PointsFile, FacesFile, CellsFile, BoundariesFile); // Set up mesh object
    std::vector<Cell> myListOfAllCls = cavity_mesh.getListOfCells(); // Get list of cells in mesh
    std::vector<Face> myListOfAllFcs = cavity_mesh.getListOfFaces(); // Get list of all faces in the mesh



    /*
        SET INITIAL DATA
    */
    double Re = 100;                                // Set Reynolds number
    double ux = 1.0, uy = 0.0;                      // Set velocity
    double L = 0.1;                                 // Set dimensions of cavity
    double nu = sqrt(ux * ux + uy * uy) * L / Re;   // Calculate nu, where, Re = u * L / nu
    int maxIterations = 150;                        // Set max number of iterations
    double alpha_P = 0.3, alpha_u = 0.7;            // under-relaxation factors



    /*
        SET UP FIELD VARIABLES AND THEIR BOUNDARY CONDITIONS 
    */
    // Velocity field
    VectorField velocity(myListOfAllCls);
    velocity.setInitialCondition(0.0, 0.0, 0.0);
    velocity.setBoundaryCondition("movingWall", "Dirichlet", 1.0, 0.0, 0.0);
    velocity.setBoundaryCondition("fixedWalls", "Dirichlet", 0.0, 0.0, 0.0);

    // Pressure field
    ScalarField pressure(myListOfAllCls);  
    pressure.setInitialCondition(0.0);
    pressure.setBoundaryCondition("movingWall", "Neumann", 0.0);
    pressure.setBoundaryCondition("fixedWalls", "Neumann", 0.0);



    /*
        FILE OUT FOR RESIDUAL NORM (need to uncomment another line in the code!)
    */
    std::ofstream outputResNorm("40x40/Results/ResNorm.txt");
    outputResNorm << "Iteration P_norm Vx_norm Vy_norm Vz_norm" << std::endl;



    /*
        FILE OUT FOR CONVERGENCE DATA
    */
    // Set up output for velocity and pressure convergence data
    std::ofstream outputConvergenceData("40x40/Results/ConvergenceData.txt");
    outputConvergenceData << "Iteration VelXAfterMomPred VelYAfterMomPred PressureAfterPressurePred VelXAfterCorr VelYAfterCorr PressureAfterCorr" << std::endl;
    Eigen::VectorXd velXAfterMomPred, velYAfterMomPred, velXAfterCorr, velYAfterCorr, presAfterPresPred, presAfterCorr;



    /* 
            SIMPLE ALGORITHM LOOP
    */
    for (int iteration = 0; iteration < maxIterations; iteration++)
    {
        /*
            SET UP AND SOLVE THE MOMENTUM EQUATION
        */
        evaluateConvectionTerm(myListOfAllFcs, velocity, velocity, "Upwind");
        evaluateDiffusionTerm(myListOfAllFcs, velocity, nu);

        // Under relax velocity
        underRelaxImplicit(velocity, alpha_u);

        // Solve the momentum equation 
        velocity.solveEquation();

        // Output velocity data after momentum predictor - for convergence data
        if (iteration > 0)
        {
            double diffX = 0.0, diffY = 0.0;
            double max_val_X = 0.0, max_val_Y = 0.0;

            for (int i = 0; i < velocity.getMatrixSize(); i++)
            {
                diffX = fabs(velocity.variableX(i) - velXAfterMomPred(i));
                diffY = fabs(velocity.variableY(i) - velYAfterMomPred(i));

                if (i == 0)
                {
                    max_val_X = diffX;
                    max_val_Y = diffY;
                }

                if (diffX > max_val_X)
                {
                    max_val_X = diffX;
                }

                if (diffY > diffY)
                {
                    max_val_Y = diffY;
                }
            }
            outputConvergenceData << iteration << " " << max_val_X << " " << max_val_Y;
        }
        velXAfterMomPred = velocity.variableX;
        velYAfterMomPred = velocity.variableY;

        // Reset diagonal of momentum coeff matrix
        for (int i = 0; i < velocity.getMatrixSize(); i++)
        {
            velocity.coeffs(i, i) = velocity.coeffs(i, i) * alpha_u;
        }



        /*
            SET UP AND SOLVE THE PRESSURE EQUATION
        */
        // Store p_old
        Eigen::VectorXd p_old = pressure.variable;  

        // Calculate F_pre and assign fluxes to faces
        for (int i = 0; i < myListOfAllFcs.size(); i++)
        {
            Face* f = &myListOfAllFcs[i];
            int idxP = f->getFaceOwner();
            int idxN = f->getFaceNeighbour();

            double F = getFlux(f, velocity, idxP, idxN);
            f->setFaceFlux(F);
        }

        // Set div u contribution to source matrix of p
        for (int i = 0; i < myListOfAllCls.size(); i++)
        {
            Cell c = myListOfAllCls[i];
            int cellIdx = c.getCellIndex();
            std::vector<Face *> myListOfFcs = c.getListOfFaces(); // Get faces for this cell
            double fluxForCell = 0.0;

            for (int j = 0; j < myListOfFcs.size(); j++)
            {
                Face* f = myListOfFcs[j]; // Pick a face from the list

                for (int m = 0; m < myListOfAllFcs.size(); m++)
                {
                    Face listFace = myListOfAllFcs[m];

                    if (f->getFaceIndex() == listFace.getFaceIndex())
                    {
                        if (listFace.getFaceOwner() == c.getCellIndex()) // face belongs to the cell
                        {
                            fluxForCell += listFace.getFaceFlux();
                        }
                        else // if the face belongs to a neighbour cell
                        {
                            fluxForCell += -listFace.getFaceFlux();
                        }
                    }
                }
            }

            pressure.source(cellIdx) += fluxForCell;
        }

        // Set up gamma for pressure field
        std::vector<double> gamma_p(myListOfAllCls.size());

        for (int i = 0; i < myListOfAllCls.size(); i++)
        {
            gamma_p[i] = 1.0 / velocity.coeffs(i, i);
        }

        // Evaluate diffusion term for pressure field and solve matrices
        evaluatePositiveDiffusionTerm(myListOfAllFcs, pressure, gamma_p); 

        // solve pressure eqn
        pressure.solveEquation();  

        // Output pressure data after pressure predictor - for convergence data
        if (iteration > 0)
        {
            double diff = 0.0, max_val = 0.0;

            for (int i = 0; i < pressure.getMatrixSize(); i++)
            {
                diff = fabs(pressure.variable(i) - presAfterPresPred(i));

                if (i == 0)
                {
                    max_val = diff;
                }

                if (diff > max_val)
                {
                    max_val = diff;
                }
            }
            outputConvergenceData << " " << max_val;
        }
        presAfterPresPred = pressure.variable;



        /*
            CORRECT FACE FLUXES
        */
        // Vector to store corrected fluxes
        std::vector<double> F_corr(myListOfAllFcs.size());

        // Calculate F_corr and assign correct face fluxes
        for (int i = 0; i < myListOfAllFcs.size(); i++)
        {
            Face f = myListOfAllFcs[i];
            int idxP = f.getFaceOwner();
            int idxN = f.getFaceNeighbour();
            double F_pre = f.getFaceFlux();
            double a_P_interpolated = 0.0;
            double fx = f.getFaceInterpolationFactor();

            // Calc corrected flux -- assumes orthogonal mesh
            if (f.getIfBoundaryFace() == true)  // If face is a boundary
            {
                F_corr[i] = F_pre; // Since boundary condition on all boundaries is grad(P) = 0 
            }
            else  // if face is not a boundary
            {
                a_P_interpolated = (fx * (1.0 / velocity.coeffs(idxP, idxP))) + ((1.0 - fx) * (1.0 / velocity.coeffs(idxN, idxN)));
                F_corr[i] = F_pre - (a_P_interpolated * f.getFaceArea() * (pressure.variable(idxN) - pressure.variable(idxP)) / f.getFaceDelta());
            }
            
            // Update old face flux
            f.setFaceFlux(F_corr[i]);
        }

        

        /*
            DIV U CHECK
        */
        //     std::vector<double> divU;
        //     for (int i = 0; i < myListOfAllCls.size(); i++)
        //     {
        //         Cell c = myListOfAllCls[i];
        //         int cellIdx = c.getCellIndex();
        //         std::vector<Face *> myListOfFcs = c.getListOfFaces(); // Get faces for this cell
        //         double fluxForCell = 0.0;

        //         for (int j = 0; j < myListOfFcs.size(); j++)
        //         {
        //             Face* f = myListOfFcs[j]; // Pick a face from the list

        //             for (int m = 0; m < myListOfAllFcs.size(); m++)
        //             {
        //                 Face listFace = myListOfAllFcs[m];

        //                 if (f->getFaceIndex() == listFace.getFaceIndex())
        //                 {
        //                     if (listFace.getFaceOwner() == c.getCellIndex()) // face belongs to the cell
        //                     {
        //                         fluxForCell += listFace.getFaceFlux();
        //                     }
        //                     else // if the face belongs to a neighbour cell
        //                     {
        //                         fluxForCell += -listFace.getFaceFlux();
        //                     }
        //                 }
        //             }
        //         }
        //         divU.push_back(fluxForCell);      
        //     }

        //     // Output final time step data to a single file
        //     std::ofstream output("40x40/Results/DivUCheck_" + std::to_string(iteration) + ".txt");

        //     output << "y_coord";
        //     for (int i = 0; i < nCellsx; i++)
        //     {
        //         output << "," << i;
        //     }
        //     output << std::endl;

        //     for (int z = 0; z < nCellsy; z++)
        //     {
        //         output << z;
        //         for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        //         {
        //             output << "," << divU[i];
        //         }
        //         output << std::endl;
        //     }
        // // END OF DIV U CHECK



        /*
            CORRECT CELL CENTRE VELOCITIES
        */
        for (int i = 0; i < myListOfAllCls.size(); i++)
        {
            Cell c = myListOfAllCls[i];
            int idx = c.getCellIndex();
            Vector grad_P = getGradientGauss(pressure, idx, myListOfAllCls);

            velocity.variableX(idx) = velocity.variableX(idx) - ((1.0 / velocity.coeffs(idx, idx)) * grad_P[0] / c.getCellVolume());
            velocity.variableY(idx) = velocity.variableY(idx) - ((1.0 / velocity.coeffs(idx, idx)) * grad_P[1] / c.getCellVolume());
            velocity.variableZ(idx) = velocity.variableZ(idx) - ((1.0 / velocity.coeffs(idx, idx)) * grad_P[2] / c.getCellVolume());
        }

        // Output velocity data after correction - for convergence data
        if (iteration > 0)
        {
            double diffX = 0.0, diffY = 0.0;
            double max_val_X = 0.0, max_val_Y = 0.0;

            for (int i = 0; i < velocity.getMatrixSize(); i++)
            {
                diffX = fabs(velocity.variableX(i) - velXAfterCorr(i));
                diffY = fabs(velocity.variableY(i) - velYAfterCorr(i));

                if (i == 0)
                {
                    max_val_X = diffX;
                    max_val_Y = diffY;
                }

                if (diffX > max_val_X)
                {
                    max_val_X = diffX;
                }

                if (diffY > diffY)
                {
                    max_val_Y = diffY;
                }
            }
            outputConvergenceData << " " << max_val_X << " " << max_val_Y;
        }
        velXAfterCorr = velocity.variableX;
        velYAfterCorr = velocity.variableY;



        /*
            UNDER-RELAX PRESSURE
        */
        underRelaxExplicit(pressure, p_old, alpha_P);

        // Output pressure data after pressure correction - for convergence data
        if (iteration > 0)
        {
            double diff = 0.0, max_val = 0.0;

            for (int i = 0; i < pressure.getMatrixSize(); i++)
            {
                diff = fabs(pressure.variable(i) - presAfterCorr(i));

                if (i == 0)
                {
                    max_val = diff;
                }

                if (diff > max_val)
                {
                    max_val = diff;
                }
            }

            outputConvergenceData << " " << max_val << std::endl;
        }
        presAfterCorr = pressure.variable;


        
        /*
            OUTPUT DATA REQUIRED FOR PLOTS FOR THIS ITERATION
        */
        outputVariableDataForPlot(nCellsx, nCellsy, "40x40/Results/velocity", velocity, iteration);
        outputVariableDataForPlot(nCellsx, nCellsy, "40x40/Results/pressure", pressure, iteration);
        outputVectorMagnitudeForPlot(nCellsx, nCellsy, "40x40/Results/velocity", velocity, iteration);
        outputNormalizedVectorCompForPlot(nCellsx, nCellsy, "40x40/Results/velocity", velocity, iteration);


        /*
            RESIDUALS 
        */
        // Calculate residuals and output to files
        calculateResidual(pressure);
        calculateResidual(velocity);
        
        // Output residual data for plots
        // outputResidualForPlot(nCellsx, nCellsy, "40x40/Results/velocity", velocity, iteration);
        // outputResidualForPlot(nCellsx, nCellsy, "40x40/Results/pressure", pressure, iteration);
        // outputAllResiduals(pressure, velocity, "40x40/Results/res", iteration);
        
        // Calculate residual norms and output to file - needed for file output (see lines 117-118)
        std::vector<double> velocityResNorm = velocity.getResidualNorm();
        double pressureResNorm = pressure.getResidualNorm();
        outputResNorm << iteration << " " << pressureResNorm << " " << velocityResNorm[0] << " " << velocityResNorm[1] << " " << velocityResNorm[2] << std::endl;



        /*
            RESET COEFF AND SOURCE MATRICES BEFORE NEXT TIME STEP
        */
        velocity.resetMatrices();
        pressure.resetMatrices();

        // Print iteration number at the end of the loop
        std::cout << "Iteration " << iteration << std::endl;
    }

    // Output at the end of simulation
    tEnd = clock();
    double simTime = double(tEnd - tStart) / double(CLOCKS_PER_SEC);
    std::cout << "Simulation complete! Time taken = " << std::setprecision(5) << simTime << " secs." << std::endl;
    
    return 0;
}