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

#include "admCustomReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admCustomReactionRate, 0);

    addToRunTimeSelectionTable
    (
        admReactionRate,
        admCustomReactionRate,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admCustomReactionRate::admCustomReactionRate
(
    admTime& runTime,
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    admReactionManager& admReacs,
    const dictionary& dict,
    const word& reactionName
)
:
    admReactionRate(admVars, admCoeffs, admReacs, dict, reactionName),

    uninhibitedValue_
    (
        admCoeffs_.addNew
        (
            dict.subDict(reactionName),
            word("rate"),
            word("custom_rate(" + reactionName + ")")
        )
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admCustomReactionRate::uninhibited(label cellIndex) const
{
    return uninhibitedValue_.evaluate(cellIndex);
}


Foam::tmp<scalarField> Foam::admCustomReactionRate::uninhibitedField() const
{
    return uninhibitedValue_.evaluateField();
}


Foam::scalar Foam::admCustomReactionRate::uninhibitedDdy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return uninhibitedValue_.ddy(wrtVar, cellIndex);
}


Foam::tmp<scalarField> Foam::admCustomReactionRate::uninhibitedDdyField
(
    const admVariable& wrtVar
) const
{
    return uninhibitedValue_.ddyField(wrtVar);
}


bool Foam::admCustomReactionRate::uninhibitedDdyNonZero
(
    const admVariable& wrtVar
) const
{
    return uninhibitedValue_.ddyNonZero(wrtVar);
}

// ************************************************************************* //
