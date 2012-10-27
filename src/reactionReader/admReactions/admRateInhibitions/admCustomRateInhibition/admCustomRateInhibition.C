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

#include "admCustomRateInhibition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admCustomRateInhibition, 0);

    addToRunTimeSelectionTable
    (
        admRateInhibition,
        admCustomRateInhibition,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admCustomRateInhibition::
    admCustomRateInhibition
(
    admTime& runTime,
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    const dictionary& dict,
    const word& inhibitionName
)
:
    admRateInhibition(admVars, admCoeffs, inhibitionName),

    inhibitionValue_
    (
        admCoeffs_.addNew
        (
            dict,
            inhibitionName,
            word("custom_" + name())
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admCustomRateInhibition::inhibition
(
    label cellIndex
) const
{
    return inhibitionValue_.evaluate(cellIndex);
}


Foam::tmp<scalarField> Foam::admCustomRateInhibition::inhibitionField() const
{
    return inhibitionValue_.evaluateField();
}


bool Foam::admCustomRateInhibition::inhibitionNonZero() const
{
    return true;
}


Foam::scalar Foam::admCustomRateInhibition::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return inhibitionValue_.ddy(wrtVar, cellIndex);
}


Foam::tmp<scalarField> Foam::admCustomRateInhibition::ddyField
(
    const admVariable& wrtVar
) const
{
    return inhibitionValue_.ddyField(wrtVar);
}


bool Foam::admCustomRateInhibition::ddyNonZero(const admVariable& wrtVar) const
{
    return inhibitionValue_.ddyNonZero(wrtVar);
}

// ************************************************************************* //
