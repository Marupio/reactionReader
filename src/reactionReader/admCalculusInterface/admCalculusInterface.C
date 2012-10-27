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

#include "admCalculusInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admCalculusInterface, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admCalculusInterface::admCalculusInterface
(
    admVariableManager& admVars,
    const word& name,
    const dimensionSet& dimensions
)
:
    runTime_(admVars.runTime()),
    mesh_(admVars.mesh()),
    admVars_(admVars),
    name_(name),
    dimensions_(dimensions)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<dimensionedScalarField>
    Foam::admCalculusInterface::evaluateDimensionedField() const
{
    return tmp<dimensionedScalarField>
    (
        new dimensionedScalarField
        (
            IOobject
            (
                evaluateName(),
                admVars_.runTime().timeName(),
                admVars_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            admVars_.mesh(),
            evaluateDims(),
            evaluateField()
        )
    );
}


Foam::tmp<scalarField> Foam::admCalculusInterface::ddyField
(
    const admVariable& wrtVar
) const
{
    return tmp<scalarField>
    (
        new scalarField(admVars_.mesh().nCells(), 0.0)
    );
}


Foam::tmp<dimensionedScalarField>
    Foam::admCalculusInterface::ddyDimensionedField
(
    const admVariable& wrtVar
) const
{
    return tmp<dimensionedScalarField>
    (
        new dimensionedScalarField
        (
            IOobject
            (
                ddyName(wrtVar),
                admVars_.runTime().timeName(),
                admVars_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            admVars_.mesh(),
            ddyDims(wrtVar),
            ddyField(wrtVar)
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
