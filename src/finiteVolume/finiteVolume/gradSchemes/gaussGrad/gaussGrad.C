/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "gaussGrad.H"
#include "extrapolatedCalculatedFvPatchField.H"




#ifdef NVTX
    #include <nvtx3/nvToolsExt.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{

    #ifdef NVTX
        nvtxRangePushA("Build Gradient");  
    #endif

        typedef typename outerProduct<vector, Type>::type GradType;
        typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

        const fvMesh& mesh = ssf.mesh();


        tmp<GradFieldType> tgGrad
        (
            new GradFieldType
            (
                IOobject
                (
                    name,
                    ssf.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<GradType>(ssf.dimensions()/dimLength, Zero),
                fvPatchFieldBase::extrapolatedCalculatedType()
            )
        );

        #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("Before First Cycle");  
    #endif


        GradFieldType& gGrad = tgGrad.ref();
        Field<GradType>& igGrad = gGrad;

        const vectorField& Sf = mesh.Sf();
        const Field<Type>& issf = ssf;



    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("First Cycle");  
    #endif

    #ifdef STDPAR

        const labelUList& owlist=mesh.lduAddr().ownerList();
        const labelUList& owstart=mesh.lduAddr().ownerStart();
        const labelUList& nelist=mesh.lduAddr().neighbourList();
        const labelUList& nestart=mesh.lduAddr().neighbourStart();

        auto iter=std::views::iota(0,igGrad.size());
        std::for_each(std::execution::par,iter.begin(),iter.end(),
                [ol=owlist.cdata(),os=owstart.cdata(),nl=nelist.cdata(),ns=nestart.cdata(),sf=Sf.cdata(),is=issf.cdata(),ig=igGrad.data()](const label& facei){
                    ig[facei]=Zero;
                    for(int i=os[facei]; i<os[facei+1];++i){
                        ig[facei]+= sf[ol[i]]*is[ol[i]];
                    }
                    for(int i=ns[facei]; i<ns[facei+1];++i){
                        ig[facei]-= sf[nl[i]]*is[nl[i]];
                    }
                });

    #else

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        forAll(owner, facei)
        {
            const GradType Sfssf = Sf[facei]*issf[facei];

            igGrad[owner[facei]] += Sfssf;
            igGrad[neighbour[facei]] -= Sfssf;
        }

    #endif

    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("Second Cycle");  
    #endif


    #ifdef STDPAR
        forAll(mesh.boundary(), patchi)
        {
        
            const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
            const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];
            
            if(mesh.boundary()[patchi].size()!=0){
            #ifdef NVTX
                nvtxRangePushA("List creation");  
            #endif
            const auto& faceIndex=mesh.boundary().facePatchIndexPatch(patchi, mesh.boundary());
            const auto& faceStart=mesh.boundary().facePatchStartPatch(patchi, mesh.boundary());

            #ifdef NVTX
                nvtxRangePop();
                nvtxRangePushA(" Cycle");  
            #endif

                std::for_each(std::execution::par,
                                std::views::iota(0).begin(),
                                std::views::iota(faceStart.size()-1).begin(), 
                                [ig=igGrad.data(),f=pFaceCells.cdata(),fir=pSf.cdata(),sec=pssf.cdata(),pstr=faceStart.cdata(),plst=faceIndex.cdata()](const label& facei){
                                    label id=f[plst[pstr[facei]]];
                                    
                                    for(int i=pstr[facei]; i<pstr[facei+1];++i){
                                        ig[id]+=fir[plst[i]]*sec[plst[i]];
                                    }
                                    
                                });

            #ifdef NVTX
                nvtxRangePop();
            #endif
                

            }
        }

    #else

        forAll(mesh.boundary(), patchi)
        {
            const labelUList& pFaceCells =
                mesh.boundary()[patchi].faceCells();

            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

            const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

            forAll(mesh.boundary()[patchi], facei)
            {
                igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
            }
        }
    
    #endif

    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("Third Cycle");  
    #endif

        igGrad /= mesh.V();

    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("Correct Boundary Conditions");  
    #endif

        gGrad.correctBoundaryConditions();

    #ifdef NVTX
        nvtxRangePop();
    #endif

    return tgGrad;
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    #ifdef NVTX
        nvtxRangePushA("calcGrad");  
    #endif

    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

        tmp<GradFieldType> tgGrad
        (
            gradf(tinterpScheme_().interpolate(vsf), name)
        );
        GradFieldType& gGrad = tgGrad.ref();

    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("Correct Boundary condition");  
    #endif


    correctBoundaryConditions(vsf, gGrad);

    #ifdef NVTX
        nvtxRangePop();
    #endif

    return tgGrad;
}


template<class Type>
void Foam::fv::gaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    auto& gGradbf = gGrad.boundaryFieldRef();

    forAll(vsf.boundaryField(), patchi)
    {
        if (!vsf.boundaryField()[patchi].coupled())
        {

        #ifdef NVTX
            nvtxRangePushA("Correct Boundary condition-1");  
        #endif
            const vectorField n
            (
                vsf.mesh().Sf().boundaryField()[patchi]
              / vsf.mesh().magSf().boundaryField()[patchi]
            );

        #ifdef NVTX
            nvtxRangePop();
            nvtxRangePushA("Correct Boundary condition-2 ");  
        #endif

            auto tmp1= vsf.boundaryField()[patchi].snGrad();

        #ifdef NVTX
            nvtxRangePop();
            nvtxRangePushA("Correct Boundary condition-3 ");  
        #endif

            gGradbf[patchi] += n *
            (
              tmp1 - (n & gGradbf[patchi])
            );

            // gGradbf[patchi] += n *
            // (
            //     vsf.boundaryField()[patchi].snGrad()
            //   - (n & gGradbf[patchi])
            // );

        #ifdef NVTX
            nvtxRangePop();
        #endif


        }
     }
}


// ************************************************************************* //
