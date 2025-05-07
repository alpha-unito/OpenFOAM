/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "lduAddressing.H"
#include "scalarField.H"


#ifdef STDPAR
    #include <atomic>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduAddressing::calcLosort() const
{
    if (losortPtr_)
    {
        FatalErrorInFunction
            << "losort already calculated"
            << abort(FatalError);
    }

    // Scan the neighbour list to find out how many times the cell
    // appears as a neighbour of the face. Done this way to avoid guessing
    // and resizing list
    labelList nNbrOfFace(size(), Foam::zero{});

    const labelUList& nbr = upperAddr();

    forAll(nbr, nbrI)
    {
        nNbrOfFace[nbr[nbrI]]++;
    }

    // Create temporary neighbour addressing
    labelListList cellNbrFaces(size());

    forAll(cellNbrFaces, celli)
    {
        cellNbrFaces[celli].setSize(nNbrOfFace[celli]);
    }

    // Reset the list of number of neighbours to zero
    nNbrOfFace = 0;

    // Scatter the neighbour faces
    forAll(nbr, nbrI)
    {
        cellNbrFaces[nbr[nbrI]][nNbrOfFace[nbr[nbrI]]] = nbrI;

        nNbrOfFace[nbr[nbrI]]++;
    }

    // Gather the neighbours into the losort array
    losortPtr_ = std::make_unique<labelList>(nbr.size(), -1);
    auto& lst = *losortPtr_;

    // Set counter for losort
    label lstI = 0;

    forAll(cellNbrFaces, celli)
    {
        const labelList& curNbr = cellNbrFaces[celli];

        forAll(curNbr, curNbrI)
        {
            lst[lstI] = curNbr[curNbrI];
            lstI++;
        }
    }
}


void Foam::lduAddressing::calcOwnerStart() const
{
    if (ownerStartPtr_)
    {
        FatalErrorInFunction
            << "owner start already calculated"
            << abort(FatalError);
    }

    const labelList& own = lowerAddr();

    ownerStartPtr_ = std::make_unique<labelList>(size() + 1, own.size());
    auto& ownStart = *ownerStartPtr_;

    // Set up first lookup by hand
    ownStart[0] = 0;
    label nOwnStart = 0;
    label i = 1;

    forAll(own, facei)
    {
        label curOwn = own[facei];

        if (curOwn > nOwnStart)
        {
            while (i <= curOwn)
            {
                ownStart[i++] = facei;
            }

            nOwnStart = curOwn;
        }
    }
}


void Foam::lduAddressing::calcLosortStart() const
{
    if (losortStartPtr_)
    {
        FatalErrorInFunction
            << "losort start already calculated"
            << abort(FatalError);
    }

    losortStartPtr_ = std::make_unique<labelList>(size() + 1, Foam::zero{});
    auto& lsrtStart = *losortStartPtr_;

    const labelList& nbr = upperAddr();

    const labelList& lsrt = losortAddr();

    // Set up first lookup by hand
    lsrtStart[0] = 0;
    label nLsrtStart = 0;
    label i = 0;

    forAll(lsrt, facei)
    {
        // Get neighbour
        const label curNbr = nbr[lsrt[facei]];

        if (curNbr > nLsrtStart)
        {
            while (i <= curNbr)
            {
                lsrtStart[i++] = facei;
            }

            nLsrtStart = curNbr;
        }
    }

    // Set up last lookup by hand
    lsrtStart[size()] = nbr.size();
}


void Foam::lduAddressing::calcLoCSR() const
{
    if (lowerCSRAddrPtr_)
    {
        FatalErrorInFunction
            << "lowerCSRAddr already calculated"
            << abort(FatalError);
    }

    lowerCSRAddrPtr_ = std::make_unique<labelList>(lowerAddr().size());
    map(lowerAddr(), *lowerCSRAddrPtr_);
}



#ifdef STDPAR


void Foam::lduAddressing::sorting_pair(labelList& index, labelList& val) const
{
    const label N=index.size();
    auto iter=std::views::iota(0,N);

    std::vector<std::pair<label,label>> tmp_pair;    
    tmp_pair.reserve(N);

    std::transform(std::execution::par_unseq,index.begin(),index.end(),val.begin(),tmp_pair.begin(),
                    [](const auto& i, const auto& v){return std::make_pair(i,v);});


    std::stable_sort(std::execution::par_unseq,tmp_pair.begin(),tmp_pair.begin()+N,
                [=](const auto& p1, const auto& p2 ){
                    return p1.second<p2.second;
                });

    
    std::for_each(std::execution::par_unseq,
                iter.begin(),
                iter.end(),
                [pr=tmp_pair.data(),id=index.data(),vl=val.data()](const auto& i){
                    id[i]=pr[i].first;
                    vl[i]=pr[i].second;
                });

}


void Foam::lduAddressing::csr_list(labelList& neigh_sort, labelList& start, const label& n, const label& N) const 
{

    label *stmp = new label[n];

    labelList ones(n+N);
    labelList nb_tmp(n+N);
    auto iter2=std::views::iota(0,n);
    
    std::fill(std::execution::par_unseq,ones.begin(),ones.end(),1);
    std::fill(std::execution::par_unseq,nb_tmp.begin(),nb_tmp.end(),0);

    std::copy(std::execution::par_unseq,neigh_sort.begin(),neigh_sort.end(),nb_tmp.begin());
    std::copy(std::execution::par_unseq,iter2.begin(),iter2.end(),nb_tmp.begin()+N);

    std::fill(std::execution::par_unseq,ones.begin()+N,ones.end(),0);
    std::fill(std::execution::par_unseq,stmp,stmp+n,0);

    sorting_pair(ones,nb_tmp);

    auto iter3=std::views::iota(0,n+N);
    std::for_each(std::execution::par,iter3.begin(),iter3.end(),
                [=,nb=nb_tmp.data(),on=ones.data()](const auto& i){

                    auto atomic_ref_stmp = std::atomic_ref<label>(stmp[nb[i]]);
                    atomic_ref_stmp.fetch_add(on[i],std::memory_order_relaxed);

                });

    std::copy(std::execution::par,stmp,stmp+n,start.begin());
    std::exclusive_scan(std::execution::par,start.begin(),start.end()+1,start.begin(),0);

    delete [] stmp;
}


void Foam::lduAddressing::sorting_pair2(labelList& index, labelList& val) 
{
    const label N=index.size();
    auto iter=std::views::iota(0,N);

    std::vector<std::pair<label,label>> tmp_pair;    
    tmp_pair.reserve(N);

    std::transform(std::execution::par_unseq,index.begin(),index.end(),val.begin(),tmp_pair.begin(),
                    [](const auto& i, const auto& v){return std::make_pair(i,v);});


    std::stable_sort(std::execution::par_unseq,tmp_pair.begin(),tmp_pair.begin()+N,
                [=](const auto& p1, const auto& p2 ){
                    return p1.second<p2.second;
                });

    
    std::for_each(std::execution::par_unseq,
                iter.begin(),
                iter.end(),
                [pr=tmp_pair.data(),id=index.data(),vl=val.data()](const auto& i){
                    id[i]=pr[i].first;
                    vl[i]=pr[i].second;
                });

}


void Foam::lduAddressing::csr_list2(const labelUList& list, labelList& lindex, labelList& start, label& n)
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
                    sorting_pair2(lindex,lsort);


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



void Foam::lduAddressing::calcownerList() const 
{
    if (ownerList_)
    {
        FatalErrorIn("lduAddressing::calcownerList() const")
            << "ownerList already calculated"
            << abort(FatalError);
    }

    const labelUList& owner=lowerAddr();
    const label N=owner.size();
    const label n=size();

    ownerList_ = std::make_unique<labelList>(N,0);
    ownerStart_ = std::make_unique<labelList>(n+1,0);

    labelList owner_sort(N);

    labelList& olist= *ownerList_;
    labelList& ostart= *ownerStart_;
    labelList& osort= owner_sort;

    std::copy(std::execution::par_unseq,std::views::iota(0).begin(),std::views::iota(N).begin(),olist.begin()); 
    std::copy(std::execution::par_unseq,owner.begin(),owner.end(),osort.begin());

    sorting_pair(olist,osort);

    csr_list(osort,ostart,n,N);

}


void Foam::lduAddressing::calcneighbourList() const
{
    if (neighbourList_)
    {
        FatalErrorIn("lduAddressing::calcneighbourList() const")
            << "neighbourList already calculated"
            << abort(FatalError);
    }

    const labelUList& neigh=upperAddr();
    const label N=neigh.size();
    const label n=size();

    neighbourList_ = std::make_unique<labelList>(N,0);
    neighbourStart_ = std::make_unique<labelList>(n+1,0);
    labelList neigh_sort(N);

    labelList& nlist= *neighbourList_;
    labelList& nstart= *neighbourStart_;
    labelList& nsort= neigh_sort;


    std::copy(std::execution::par_unseq,std::views::iota(0).begin(),std::views::iota(N).begin(),nlist.begin());

    std::copy(std::execution::par_unseq,neigh.begin(),neigh.end(),nsort.begin());
    sorting_pair(nlist,nsort);
    
    csr_list(nsort,nstart,n,N);

}




#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::lduAddressing::losortAddr() const
{
    if (!losortPtr_)
    {
        calcLosort();
    }

    return *losortPtr_;
}


const Foam::labelUList& Foam::lduAddressing::ownerStartAddr() const
{
    if (!ownerStartPtr_)
    {
        calcOwnerStart();
    }

    return *ownerStartPtr_;
}


const Foam::labelUList& Foam::lduAddressing::losortStartAddr() const
{
    if (!losortStartPtr_)
    {
        calcLosortStart();
    }

    return *losortStartPtr_;
}


const Foam::labelUList& Foam::lduAddressing::lowerCSRAddr() const
{
    if (!lowerCSRAddrPtr_)
    {
        calcLoCSR();
    }

    return *lowerCSRAddrPtr_;
}


#ifdef STDPAR

const Foam::labelUList& Foam::lduAddressing::ownerList() const
{
    if (!ownerList_)
    {
        calcownerList();
    }

    return *ownerList_;
}

const Foam::labelUList& Foam::lduAddressing::ownerStart() const 
{
    if (!ownerStart_)
    {
        calcownerList();
    }

    return *ownerStart_;
}

const Foam::labelUList& Foam::lduAddressing::neighbourList() const
{
    if (!neighbourList_)
    {
        calcneighbourList();
    }

    return *neighbourList_;
}

const Foam::labelUList& Foam::lduAddressing::neighbourStart() const
{
    if (!neighbourStart_)
    {
        calcneighbourList();
    }

    return *neighbourStart_;
}

#endif

#ifdef STDPAR
void Foam::lduAddressing::clearOut()
{
    losortPtr_.reset(nullptr);
    ownerStartPtr_.reset(nullptr);
    losortStartPtr_.reset(nullptr);
    lowerCSRAddrPtr_.reset(nullptr);

    ownerList_.reset(nullptr);
    ownerStart_.reset(nullptr);
    neighbourList_.reset(nullptr);
    neighbourStart_.reset(nullptr);
}
#else
void Foam::lduAddressing::clearOut()
{
    losortPtr_.reset(nullptr);
    ownerStartPtr_.reset(nullptr);
    losortStartPtr_.reset(nullptr);
    lowerCSRAddrPtr_.reset(nullptr);
}
#endif

Foam::label Foam::lduAddressing::triIndex(const label a, const label b) const
{
    label own = min(a, b);

    label nbr = max(a, b);

    label startLabel = ownerStartAddr()[own];

    label endLabel = ownerStartAddr()[own + 1];

    const labelUList& neighbour = upperAddr();

    for (label i=startLabel; i<endLabel; i++)
    {
        if (neighbour[i] == nbr)
        {
            return i;
        }
    }

    // If neighbour has not been found, something has gone seriously
    // wrong with the addressing mechanism
    FatalErrorInFunction
        << "neighbour " << nbr << " not found for owner " << own << ". "
        << "Problem with addressing"
        << abort(FatalError);

    return -1;
}


Foam::Tuple2<Foam::label, Foam::scalar> Foam::lduAddressing::band() const
{
    const labelUList& owner = lowerAddr();
    const labelUList& neighbour = upperAddr();

    labelList cellBandwidth(size(), Foam::zero{});

    forAll(neighbour, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // Note: mag not necessary for correct (upper-triangular) ordering.
        label diff = nei-own;
        cellBandwidth[nei] = max(cellBandwidth[nei], diff);
    }

    label bandwidth = max(cellBandwidth);

    // Do not use field algebra because of conversion label to scalar
    scalar profile = 0.0;
    forAll(cellBandwidth, celli)
    {
        profile += 1.0*cellBandwidth[celli];
    }

    return Tuple2<label, scalar>(bandwidth, profile);
}


// ************************************************************************* //
