/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "PtrListOps.H"

#ifdef STDPAR

    #include<atomic>

#endif

#ifdef STDPAR

void Foam::fvBoundaryMesh::sorting_pair(labelList& index, labelList& val) const
{
    const label N=index.size();
    auto iter=std::views::iota(0,N);

    std::vector<std::pair<label,label>> tmp_pair;    
    tmp_pair.reserve(N);

    std::transform(std::execution::par,index.begin(),index.end(),val.begin(),tmp_pair.begin(),
                    [](const auto& i, const auto& v){return std::make_pair(i,v);});


    std::stable_sort(std::execution::par,tmp_pair.begin(),tmp_pair.begin()+N,
                [=](const auto& p1, const auto& p2 ){
                    return p1.second<p2.second;
                });

    
    std::for_each(std::execution::par,
                iter.begin(),
                iter.end(),
                [pr=tmp_pair.data(),id=index.data(),vl=val.data()](const auto& i){
                    id[i]=pr[i].first;
                    vl[i]=pr[i].second;
                });

}


void Foam::fvBoundaryMesh::csr_list2(const labelUList& list, labelList& lindex, labelList& start, label& n)const
{ 

    const label N=list.size();
    labelList lsort;
    labelList ones;

    lindex.resize(N);
    lsort.resize(N);
    ones.resize(N);

    auto iter=std::views::iota(0,N);


    std::copy(std::execution::par_unseq,list.begin(),list.end(),lsort.begin());
    std::copy(std::execution::par_unseq,iter.begin(),iter.end(),lindex.begin());
    std::fill(std::execution::par_unseq,ones.begin(),ones.end(),1);
    sorting_pair(lindex,lsort);


    std::for_each(std::execution::par,iter.begin(),iter.end()-1,
            [ns=lsort.data(),ts=ones.data()](const auto& x){
                if(ns[x]==ns[x+1]){
                    ts[x+1]=0;
                }
            });
    std::inclusive_scan(std::execution::par_unseq,ones.begin(),ones.end(),ones.begin());

    n=ones.last();
    start.resize(n+1);

    std::fill(std::execution::par, start.data(), start.data()+n+1, 0);
    std::for_each(std::execution::par, iter.begin(), iter.end(),
        [s=start.data(), ts = ones.cdata()](const auto& x) {
            std::atomic_ref<label>(s[ts[x] - 1]).fetch_add(1, std::memory_order_relaxed);
        });

    std::exclusive_scan(std::execution::par_unseq,start.begin(),start.end(),start.begin(),0);

}


void Foam::fvBoundaryMesh::calcfacePatchStartANDIndex(const Foam::fvBoundaryMesh& bound) const {

    // nbound_=bound.size();
    facePatchIndex_= new List<labelList>(bound.size());
    facePatchStart_= new List<labelList>(bound.size());

    auto& faceIndexRef=*facePatchIndex_;
    auto& faceStartRef=*facePatchStart_;

    forAll(bound, patchi)
        {
        
            const labelUList& pFaceCells = bound[patchi].faceCells();


            if(bound[patchi].size()!=0){
            
                label n=0;

                csr_list2(pFaceCells,faceIndexRef[patchi],faceStartRef[patchi],n);

                // labelList faceIndex;
                // labelList faceStart;
                // label n=0;
                // csr_list2(pFaceCells,faceIndex,faceStart,n);

                // for(int facei=0; facei<n+1; ++facei){
                //     faceStartRef[patchi].append(faceStart[facei]);
                // }

                // for(int facei=0; facei<faceIndex.size(); ++facei){
                //     faceIndexRef[patchi].append(faceIndex[facei]);	
                // }


            }

        }

}

#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::addPatches(const polyBoundaryMesh& pbm)
{
    // Set boundary patches
    fvPatchList& patches = *this;

    patches.resize_null(pbm.size());

    forAll(patches, patchi)
    {
        patches.set(patchi, fvPatch::New(pbm[patchi], *this));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m
)
:
    fvPatchList(),
    mesh_(m)
{}


Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m,
    const polyBoundaryMesh& pbm
)
:
    fvPatchList(),
    mesh_(m)
{
    addPatches(pbm);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvBoundaryMesh::indices
(
    const wordRe& matcher,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(matcher, useGroups);
}


Foam::labelList Foam::fvBoundaryMesh::indices
(
    const wordRes& matcher,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(matcher, useGroups);
}


Foam::labelList Foam::fvBoundaryMesh::indices
(
    const wordRes& select,
    const wordRes& ignore,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(select, ignore, useGroups);
}


Foam::label Foam::fvBoundaryMesh::findPatchID(const word& patchName) const
{
    if (patchName.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, patchName);
}


void Foam::fvBoundaryMesh::movePoints()
{
    fvPatchList& patches = *this;

    for (fvPatch& p : patches)
    {
        p.initMovePoints();
    }

    for (fvPatch& p : patches)
    {
        p.movePoints();
    }
}


Foam::UPtrList<const Foam::labelUList>
Foam::fvBoundaryMesh::faceCells() const
{
    const fvPatchList& patches = *this;

    UPtrList<const labelUList> list(patches.size());

    forAll(list, patchi)
    {
        list.set(patchi, &patches[patchi].faceCells());
    }

    return list;
}


Foam::lduInterfacePtrsList Foam::fvBoundaryMesh::interfaces() const
{
    const fvPatchList& patches = *this;

    lduInterfacePtrsList list(patches.size());

    forAll(list, patchi)
    {
        const lduInterface* lduPtr = isA<lduInterface>(patches[patchi]);

        if (lduPtr)
        {
            list.set(patchi, lduPtr);
        }
    }

    return list;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::readUpdate(const polyBoundaryMesh& pbm)
{
    addPatches(pbm);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::fvPatch& Foam::fvBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


Foam::fvPatch& Foam::fvBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


// ************************************************************************* //
