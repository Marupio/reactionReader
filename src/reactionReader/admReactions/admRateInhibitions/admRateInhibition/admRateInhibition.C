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

#include "admRateInhibition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admRateInhibition, 0);

    defineRunTimeSelectionTable(admRateInhibition, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admRateInhibition::admRateInhibition
(
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    const word& name
)
:
    admCalculusInterface
    (
        admVars,
        name,
        dimless
    ),

    admCoeffs_(admCoeffs),
    
    lastValue_(mesh_.nCells()),
    lastMilestone_(-1)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admRateInhibition::evaluate(label cellIndex) const
{
    if (lastMilestone_ != runTime_.milestone())
    {
        return inhibition(cellIndex);
    }
    return lastValue_[cellIndex];
}


Foam::tmp<scalarField> Foam::admRateInhibition::evaluateField() const
{
    if (lastMilestone_ != runTime_.milestone())
    {
        lastValue_ = inhibitionField();
        lastMilestone_ = runTime_.milestone();
    }
    return tmp<scalarField>
    (
        new scalarField(lastValue_)
    );
}


bool Foam::admRateInhibition::evaluateNonZero() const
{
    return inhibitionNonZero();
}


bool Foam::admRateInhibition::inhibitionNonZero() const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newAdmRateInhibition.C"

// ************************************************************************* //
