/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

Description
    lduMatrix member operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::sumDiag()
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
    scalarField& Diag = diag();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face=0; face<l.size(); face++)
    {
        Diag[l[face]] += Lower[face];
        Diag[u[face]] += Upper[face];
    }
}


void Foam::lduMatrix::negSumDiag()
{

    #ifdef STDPAR

        scalarField& Diag = diag();
        const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
        const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();

        const labelUList& owlist=lduAddr().ownerList();
        const labelUList& owstart=lduAddr().ownerStart();
        const labelUList& nelist=lduAddr().neighbourList();
        const labelUList& nestart=lduAddr().neighbourStart();

        auto iter=std::views::iota(0,Diag.size());
        std::for_each(std::execution::par,iter.begin(),iter.end(),
                [ol=owlist.cdata(),os=owstart.cdata(),nl=nelist.cdata(),ns=nestart.cdata(),Lw=Lower.cdata(),Up=Upper.cdata(),D=Diag.data()](const label& facei){
                    for(int i=os[facei]; i<os[facei+1];++i){
                        D[facei]-=Lw[ol[i]];
                    }
                    for(int i=ns[facei]; i<ns[facei+1];++i){
                        D[facei]-=Up[nl[i]];
                    }
                });

    #else

        const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
        const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
        scalarField& Diag = diag();

        const labelUList& l = lduAddr().lowerAddr();
        const labelUList& u = lduAddr().upperAddr();

        for (label face=0; face<l.size(); face++)
        {
            Diag[l[face]] -= Lower[face];
            Diag[u[face]] -= Upper[face];
        }

    #endif


}


void Foam::lduMatrix::sumMagOffDiag
(
    scalarField& sumOff
) const
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face = 0; face < l.size(); face++)
    {
        sumOff[u[face]] += mag(Lower[face]);
        sumOff[l[face]] += mag(Upper[face]);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::lduMatrix::operator=(const lduMatrix& A)
{
    if (this == &A)
    {
        return;  // Self-assignment is a no-op
    }

    if (A.hasLower())
    {
        lower() = A.lower();
    }
    else
    {
        lowerPtr_.reset(nullptr);
    }

    if (A.hasUpper())
    {
        upper() = A.upper();
    }
    else
    {
        upperPtr_.reset(nullptr);
    }

    if (A.hasDiag())
    {
        diag() = A.diag();
    }
}


void Foam::lduMatrix::operator=(lduMatrix&& A)
{
    if (this == &A)
    {
        return;  // Self-assignment is a no-op
    }

    diagPtr_ = std::move(A.diagPtr_);
    upperPtr_ = std::move(A.upperPtr_);
    lowerPtr_ = std::move(A.lowerPtr_);
}


void Foam::lduMatrix::negate()
{
    if (diagPtr_)
    {
        diagPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }
}


void Foam::lduMatrix::operator+=(const lduMatrix& A)
{
    if (A.hasDiag())
    {
        diag() += A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        upper() += A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() += A.upper();
        lower() += A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.hasUpper())
        {
            lower() += A.upper();
            upper() += A.upper();
        }
        else
        {
            lower() += A.lower();
            upper() += A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() += A.lower();
        upper() += A.upper();
    }
    else if (diagonal())
    {
        if (A.hasUpper())
        {
            upper() = A.upper();
        }

        if (A.hasLower())
        {
            lower() = A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        if (debug > 1)
        {
            WarningInFunction
                << "Unknown matrix type combination" << nl
                << "    this : " << this->matrixTypeName()
                << "    A    : " << A.matrixTypeName() << endl;
        }
    }
}


void Foam::lduMatrix::operator-=(const lduMatrix& A)
{

    #ifdef NVTX
        nvtxRangePushA("Phase 1");  
    #endif
    
    if (A.diagPtr_)
    {
        diag() -= A.diag();
    }

    #ifdef NVTX
        nvtxRangePop();
        nvtxRangePushA("Phase 2");  
    #endif
    
    if (symmetric() && A.symmetric())
    {
        upper() -= A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() -= A.upper();
        lower() -= A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.hasUpper())
        {
            lower() -= A.upper();
            upper() -= A.upper();
        }
        else
        {
            lower() -= A.lower();
            upper() -= A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() -= A.lower();
        upper() -= A.upper();
    }
    else if (diagonal())
    {
        if (A.hasUpper())
        {

        #ifdef STDPAR

            std::transform(std::execution::par_unseq,A.upper().begin(),A.upper().end(),upper().begin(),[](const auto& elem) { return -elem; });
            
        #else
            upper() = -A.upper();

        #endif

        }

        if (A.hasLower())
        {
            lower() = -A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        if (debug > 1)
        {
            WarningInFunction
                << "Unknown matrix type combination" << nl
                << "    this : " << this->matrixTypeName()
                << "    A    : " << A.matrixTypeName() << endl;
        }
    }

    #ifdef NVTX
        nvtxRangePop();
    #endif

}


void Foam::lduMatrix::operator*=(const scalarField& sf)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    // Non-uniform scaling causes a symmetric matrix
    // to become asymmetric
    if (symmetric() || asymmetric())
    {
        scalarField& upper = this->upper();
        scalarField& lower = this->lower();

        const labelUList& l = lduAddr().lowerAddr();
        const labelUList& u = lduAddr().upperAddr();

        for (label face=0; face<upper.size(); face++)
        {
            upper[face] *= sf[l[face]];
        }

        for (label face=0; face<lower.size(); face++)
        {
            lower[face] *= sf[u[face]];
        }
    }
}


void Foam::lduMatrix::operator*=(scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }
}


// ************************************************************************* //
