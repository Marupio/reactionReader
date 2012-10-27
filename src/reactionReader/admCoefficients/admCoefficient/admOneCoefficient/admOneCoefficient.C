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

#include "admOneCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admOneCoefficient, 0);

    addToRunTimeSelectionTable
    (
        admCoefficient,
        admOneCoefficient,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admOneCoefficient::admOneCoefficient
(
    admVariableManager& admVars,
    const word& coeffName,
    const dimensionSet& dimensions
)
:
    admCoefficient(admVars, coeffName, dimless)
{}


Foam::admOneCoefficient::admOneCoefficient
(
    admTime& runTime,
    const fvMesh& mesh,
    admVariableManager& admVars,
    const dictionary& dict,
    const word& coeffName
)
:
    admCoefficient(admVars, coeffName, dimless)
{
    if (dict.subDict(coeffName).found("dimensions"))
    {
        dimensions_.reset(dict.subDict(coeffName).lookup("dimensions"));
    }
}


// ************************************************************************* //
