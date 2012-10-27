/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "admImplicitVariable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const dimensionedScalar& admImplicitVariable::gamma() const
{
    return gamma_;
}


inline dimensionedScalar& admImplicitVariable::gamma()
{
    return gamma_;
}


inline const scalar& admImplicitVariable::convergence() const
{
    return convergence_;
}


inline const scalar& admImplicitVariable::minErrorScale() const
{
    return minErrorScale_;
}


inline const Switch& admImplicitVariable::autoSolve() const
{
    return autoSolve_;
}


inline scalar& admImplicitVariable::autoSolveConvergence()
{
    return autoSolveConvergence_;
}


inline const scalar& admImplicitVariable::autoSolveConvergence() const
{
    return autoSolveConvergence_;
}


inline const label& admImplicitVariable::autoSolveMaxIter() const
{
    return autoSolveMaxIter_;
}


inline const label& admImplicitVariable::autoSolveIndex() const
{
    return autoSolveIndex_;
}


inline label& admImplicitVariable::autoSolveIndex()
{
    return autoSolveIndex_;
}


const volScalarField& admImplicitVariable::operator()() const
{
    return dataPtr_();
}


volScalarField& admImplicitVariable::operator()()
{
    return dataPtr_();
}


inline const Foam::scalar& Foam::admImplicitVariable::operator()
(
    label cellIndex
) const
{
    return dataPtr_().internalField()[cellIndex];
}


inline Foam::scalar& Foam::admImplicitVariable::operator()(label cellIndex)
{
    return dataPtr_().internalField()[cellIndex];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
