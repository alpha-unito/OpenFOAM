/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "extrapolatedCalculatedFvPatchFields.H"



#ifdef NVTX
    #include <nvtx3/nvToolsExt.h>
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrate
(
    Field<Type>& ivf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{


    #ifdef NVTX
        nvtxRangePushA("first");  
    #endif

    const fvMesh& mesh = ssf.mesh();

    const Field<Type>& issf = ssf;

    #ifdef STDPAR

        const labelUList& owlist=mesh.lduAddr().ownerList();
        const labelUList& owstart=mesh.lduAddr().ownerStart();
        const labelUList& nelist=mesh.lduAddr().neighbourList();
        const labelUList& nestart=mesh.lduAddr().neighbourStart();

        auto iter=std::views::iota(0,ivf.size());
        std::for_each(std::execution::par,iter.begin(),iter.end(),
                [ol=owlist.cdata(),os=owstart.cdata(),nl=nelist.cdata(),ns=nestart.cdata(),is=issf.cdata(),ig=ivf.data()](const label& facei){
                    for(int i=os[facei]; i<os[facei+1];++i){
                        ig[facei]+= is[ol[i]];
                    }
                    for(int i=ns[facei]; i<ns[facei+1];++i){
                        ig[facei]-= is[nl[i]];
                    }
                }); 

    #else

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        forAll(owner, facei)
        {
            ivf[owner[facei]] += issf[facei];
            ivf[neighbour[facei]] -= issf[facei];
        }

    #endif


    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("second");  
    #endif

#ifdef STDPAR

    forAll(mesh.boundary(), patchi)
    {

        if(mesh.boundary()[patchi].size()!=0){

        const labelUList& pFaceCells =mesh.boundary()[patchi].faceCells();
        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        const auto& faceIndex=mesh.boundary().facePatchIndexPatch(patchi, mesh.boundary());
        const auto& faceStart=mesh.boundary().facePatchStartPatch(patchi, mesh.boundary());

        auto iter=std::views::iota(0,faceStart.size()-1);
        std::for_each(std::execution::par_unseq,
                        iter.begin(),
                        iter.end(), 
                        [iv=ivf.data(),pf=pFaceCells.cdata(),sec=pssf.cdata(),pstr=faceStart.cdata(),plst=faceIndex.cdata()](const label& facei){
                            label id=pf[plst[pstr[facei]]];
                            for(int i=pstr[facei]; i<pstr[facei+1];++i){
                                iv[id]+=sec[plst[i]];
                            }
                        });
        }

    }

#else

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
    }

#endif


    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("third");  
    #endif

    ivf /= mesh.Vsc();


    #ifdef NVTX
        nvtxRangePop();
    #endif


}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{

    #ifdef NVTX
        nvtxRangePushA("build surfaceIntegrate");  
    #endif

    const fvMesh& mesh = ssf.mesh();


    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions()/dimVol, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );

    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("surf");  
    #endif


    surfaceIntegrate(vf.primitiveFieldRef(), ssf);

    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("coorect");  
    #endif

    vf.correctBoundaryConditions();


    #ifdef NVTX
        nvtxRangePop();
    #endif

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        fvc::surfaceIntegrate(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions(), Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += ssf[facei];
        vf[neighbour[facei]] += ssf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vf[pFaceCells[facei]] += pssf[facei];
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> surfaceSum
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf = surfaceSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
