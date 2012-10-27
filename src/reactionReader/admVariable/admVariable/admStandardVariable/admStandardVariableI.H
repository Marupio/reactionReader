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

#include "admStandardVariable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const dimensionedScalar& admStandardVariable::gamma() const
{
    return gamma_;
}


inline dimensionedScalar& admStandardVariable::gamma()
{
    return gamma_;
}


inline const scalar& admStandardVariable::convergence() const
{
    return convergence_;
}


inline const scalar& admStandardVariable::minErrorScale() const
{
    return minErrorScale_;
}


inline const Switch& admStandardVariable::changedByUdf() const
{
    return changedByUdf_;
}


const volScalarField& admStandardVariable::operator()() const
{
    return dataPtr_();
}


volScalarField& admStandardVariable::operator()()
{
    return dataPtr_();
}


inline const Foam::scalar& Foam::admStandardVariable::operator()
(
    label cellIndex
) const
{
    return dataPtr_().internalField()[cellIndex];
}


inline Foam::scalar& Foam::admStandardVariable::operator()(label cellIndex)
{
    return dataPtr_().internalField()[cellIndex];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
