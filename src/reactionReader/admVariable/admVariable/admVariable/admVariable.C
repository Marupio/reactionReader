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

#include "admVariable.H"
#include "admVariableManager.H"
#include "IOstreams.H"
#include "token.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admVariable, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admVariable::admVariable
(
    admTime& runTime,
    admVariableManager& admVars,
    const dictionary& dict,
    const label globalIndex,
    const label localIndex
)
:
    dictionary(dict),
    runTime_(runTime),
    mesh_(admVars.mesh()),
    admVars_(admVars),
    globalIndex_(globalIndex),
    localIndex_(localIndex),
    varType_(vtnone)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::admVariable::varTypeEnum Foam::admVariable::wordToVarType
(
    const Foam::word& wordIn
)
{
    if (wordIn == "standard")
    {
        return vtstandard;
    }
    else if (wordIn == "implicit")
    {
        return vtimplicit;
    }
    else if (wordIn == "derived")
    {
        return vtderived;
    }
    else // wordIn == "none" or incorrect word
    {
        return vtnone;
    }
}


Foam::word Foam::admVariable::varTypeToWord
(
    const Foam::admVariable::varTypeEnum& vtIn
)
{
    switch (vtIn)
    {
        case vtstandard:
            return "standard";
        case vtimplicit:
            return "implicit";
        case vtderived:
            return "derived";
        case vtnone:
            return "none";
    }
    FatalErrorIn("admVariable::varTypeToWord")
        << "Unknown variable type (is " << vtIn << ", should be 0-3)"
        << abort(FatalError);
    return "error";
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::admVariable::operator==(const Foam::admVariable& as) const
{
    // Only check global index & trust the object is the same
    return (globalIndex_ == as.globalIndex_);
    // TODO If the design requires only one copy of each variable can exist,
    // then I should use a static pointer list, and == checks if it is itself
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
