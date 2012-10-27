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

#include "admConstantCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admConstantCoefficient, 0);

    addToRunTimeSelectionTable
    (
        admCoefficient,
        admConstantCoefficient,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admConstantCoefficient::admConstantCoefficient
(
    admTime& runTime,
    const fvMesh& mesh,
    admVariableManager& admVars,
    const dictionary& dict,
    const word& coeffName
)
:
    admCoefficient(admVars, coeffName, dimless),
    value_(0)
{
    // Two configurations to check for:
    //  Shorthand configuration:
    //      coeffName   [entry];
    //  Full configuration:
    //      coeffName
    //      {
    //          type    constant;
    //          value   [entry];
    //      }
    //  [entry] can have 6 formats:
    //      ignoredWord [dimensionSet] scalar;
    //      [dimensionSet] scalar;
    //      scalar;
    //      ignoredWord [dimensionSet] "equation";
    //      [dimensionSet] "equation";
    //      "equation";
    //
    //  "equation" formats can not vary spatially or temporally
    word lookupName(coeffName);
    dictionary const * dictPtr(&dict);

    if (dict.isDict(coeffName))
    {
        // Full configuration (sub-dictionary detected)
        lookupName = "value";
        dictPtr = &dict.subDict(coeffName);
    }

    dimensionedScalar ds
    (
        runTime_.readDimensionedScalarOrEquation
        (
            dictPtr->lookup(lookupName),
            coeffName
        )
    );
    value_ = ds.value();
    dimensions_.reset(ds.dimensions());
/*
    if (equationReader::isDimensionedScalar(dictPtr->lookup(lookupName)))
    {
        dimensionedScalar ds(dictPtr->lookup(lookupName));
        value_ = ds.value();
        dimensions_.reset(ds.dimensions());
    }
    else if
    (
        equationReader::isNamelessDimensionedScalar
        (
            dictPtr->lookup(lookupName)
        )
    )
    {
        ITstream is(dictPtr->lookup(lookupName));
        dimensions_.reset(dimensionSet(is));
        value_ = readScalar(is);
    }
    else if (equationReader::isScalar(dictPtr->lookup(lookupName)))
    {
        value_ = readScalar(is);
    }
    else if (equationReader::isEquation(dictPtr->lookup(lookupName)))
    {
        word eqnName("_once_" + coeffName);
        equation eqn(is, eqnName);
        equationReader& eqns(runTime_.eqns());
        eqns.readEquation(eqn);
        dimensionedScalar ds(eqns.evaluateDimensioned(eqnName));
        value_ = ds.value();
        dimensions_.reset(ds.dimensions());
    }
    else
    {
        FatalErrorIn("admConstantCoefficient::admConstantCoefficient")
            << "Coefficient " << coeffName << " dictionary entry is invalid. "
            << "Valid formats are:" << token::NL << token::TAB
            << "ignoredWord [dimensionSet] scalar;" << token::NL << token::TAB
            << "[dimensionSet] scalar;" << token::NL << token::TAB
            << "scalar;" << token::NL << token::TAB
            << "ignoredWord [dimensionSet] \"equation\";" << token::NL
            << token::TAB << "[dimensionSet] \"equation\";" << token::NL
            << token::TAB << "\"equation\";"
            << abort(FatalError);
    }
*/
}

// ************************************************************************* //
