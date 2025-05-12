/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "\nCalculating temperature distribution\n" << endl;

    scalar MatrixConstructionTime=0.0;
    scalar SolverTime=0.0;
    scalar timefirstAssembly=0.0;
    lduMatrix::debug=0;
    label iter=0;
    scalar startTime = runTime.elapsedCpuTime();

    while (simple.loop())
    {
        while (simple.correctNonOrthogonal())
        {

            #ifdef NVTX
                nvtxRangePushA("Matrix construction");  
            #endif

            scalar time1=runTime.elapsedCpuTime();
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );
            MatrixConstructionTime += runTime.elapsedCpuTime() - time1;

            scalar time2=runTime.elapsedCpuTime();
            #ifdef NVTX
                nvtxRangePop();
                nvtxRangePushA("Constraints");  
            #endif
            fvOptions.constrain(TEqn);
            #ifdef NVTX
                nvtxRangePop();
                nvtxRangePushA("Solver");  
            #endif
            TEqn.solve();
            #ifdef NVTX
                nvtxRangePop();
                nvtxRangePushA("Correct");  
            #endif
            fvOptions.correct(T);
            #ifdef NVTX
                nvtxRangePop();
            #endif
            SolverTime += runTime.elapsedCpuTime() - time2;

        }

        if(iter==0){
            timefirstAssembly=MatrixConstructionTime;
        }
        ++iter;
        
        // #ifdef NVTX
        //     nvtxRangePushA("Write Fields");  
        //     #include "write.H"
        //     runTime.printExecutionTime(Info);
        //     nvtxRangePop();
        // #endif

    }


    scalar elapsed = runTime.elapsedCpuTime() - startTime;
    Info << "Matrix construction first time: " << timefirstAssembly << " s" << endl;
    Info << "Matrix construction time: " << MatrixConstructionTime << " s" << endl;
    Info << "Solver time: " << SolverTime << " s" << endl;
    Info << "Total Cycle time: " << elapsed << " s" << endl;

    runTime.printExecutionTime(Info);
    Info<< "End\n" << endl;

    return 0;


}


// ************************************************************************* //
