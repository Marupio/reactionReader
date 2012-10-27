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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::admCoefficient> Foam::admCoefficient::New
(
    admTime& runTime,
    const fvMesh& mesh,
    admVariableManager& admVars,
    const dictionary& dict,
    const word& coeffName
)
{
    // Check for shorthand entry:
    //      coeffName   value;
    // (defaults to "constant" type)
    // Full format uses a sub-dictionary:
    //      coeffName
    //      {
    //          entries     values; ... etc.
    //      }

    word coeffType("constant");

    if (dict.isDict(coeffName))
    {
        // Full format (sub-dictionary detected)
        coeffType = word(dict.subDict(coeffName).lookup("type"));
    }

    if (debug)
    {
        Info<< "admCoefficient::New(): type = " << coeffType
            << ", name = " << coeffName << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(coeffType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "admCoefficient::New"
        )   << "Coefficient " << coeffName << " has an unknown coefficient "
            << "type " << coeffType << endl << endl
            << "Valid coefficient types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return cstrIter()(runTime, mesh, admVars, dict, coeffName);
}


// ************************************************************************* //
