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

#include "admCoefficientManager.H"
#include "IOstreams.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admCoefficientManager, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admCoefficientManager::readDict()
{
    wordList coeffNames(admCoefficientDict_.toc());
    
    forAll(coeffNames, i)
    {
        autoPtr<admCoefficient> tAdmCoeff
        (
            admCoefficient::New
            (
                runTime_,
                mesh_,
                admVars_,
                admCoefficientDict_,
                coeffNames[i]
            )
        );

        insert
        (
            coeffNames[i],
            tAdmCoeff.ptr()
        );
    }

    // Add coefficients to equationReader data sources
    for
    (
        HashPtrTable<admCoefficient>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        runTime_.eqns().addSource(* iter());
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admCoefficientManager::admCoefficientManager
(
    admTime& runTime,
    const fvMesh& mesh,
    admVariableManager& admVars,
    const word& coefficientDictName
)
:
    HashPtrTable<admCoefficient>(1000),

    runTime_(runTime),
    
    mesh_(mesh),
    
    admVars_(admVars),
    
    admCoefficientDict_
    (
        IOobject
        (
            coefficientDictName,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    zero_
    (
        admVars,
        "zero",
        dimless
    ),
    
    one_
    (
        admVars,
        "one",
        dimless
    )
{
    readDict();
    admVars_.initialize();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::admCoefficientManager::toc() const
{
    return HashPtrTable<admCoefficient>::toc();
/*    wordList returnMe(0);
    for
    (
        HashPtrTable<admCoefficient>::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        returnMe.setSize(returnMe.size() + 1);
        returnMe[returnMe.size() - 1] = iter()->name();
    }
    return returnMe;*/
}


const Foam::admCoefficient& Foam::admCoefficientManager::operator()
(
    const word& coeffName
) const
{
    admCoefficient * coeffPtr(this->operator[](coeffName));
#   ifdef FULLDEBUG
    if (coeffPtr == NULL)
    {
        FatalErrorIn("admCoefficientManager::operator()")
            << coeffName << " is not the name of a coefficient.  Valid names "
            << "are: " << toc()
            << abort(FatalError);
    }
#   endif
    return * coeffPtr;
}


const admCoefficient& Foam::admCoefficientManager::addNew
(
    const dictionary& dict,
    const word& lookupName,
    const word& newCoeffName
)
{
    word coeffName(newCoeffName);
    if (newCoeffName == word::null)
    {
        coeffName = lookupName;
    }
    if (found(coeffName))
    {
        FatalErrorIn("admCoefficientManager::addNew")
            << "Attempting to add duplicate " << coeffName << " coefficient."
            << abort(FatalError);
    }
    dictionary tempDict(dict);
    tempDict.changeKeyword(lookupName, coeffName);
    autoPtr<admCoefficient> tAdmCoeff
    (
        admCoefficient::New
        (
            runTime_,
            mesh_,
            admVars_,
            tempDict,
            coeffName
        )
    );

    insert
    (
        coeffName,
        tAdmCoeff.ptr()
    );
    runTime_.eqns().addSource
    (
        operator()(coeffName)
    );
    
    return operator()(coeffName);
}


const admCoefficient& Foam::admCoefficientManager::readCoefficient
(
    const dictionary& dict,
    const word& lookupName,
    const word& newCoeffName
)
{
    word coeffName(newCoeffName);
    if (newCoeffName == word::null)
    {
        coeffName = lookupName;
    }

    // Check dictionary format - can be any of:
    //      lookupName      {subDict}  (full inline coefficient definition)
    //      lookupName      word; (gives coeffName)
    //      lookupName      dimensionedScalarOrEqutaion; (short for constant
    //                      coefficient)

    const admCoefficient * coeffPtr;

    if
    (
        dict.isDict(lookupName)
     || !equationReader::isWord
        (
            dict.lookup(lookupName)
        )
    )
    {
        // Either a full inline description (subdict), or a short format for
        // a constant coefficient.  Each requires a new coefficient to be
        // created.
    
        coeffPtr = & addNew
        (
            dict,
            lookupName,
            coeffName
        );

        /*dictionary tempDict;
        tempDict.set(newCoeffName, dict.subDict(lookupName));
        coeffPtr = & addNew
        (
            tempDict,
            coeffName
        );*/
    }
    else if
    (
        equationReader::isWord
        (
            dict.lookup(lookupName)
        )
    )
    {
        // Names a coefficient
        word coeffName(dict.lookup(lookupName));
        coeffPtr = & operator()(word(dict.lookup(lookupName)));
    }
    else
    {
        FatalIOErrorIn("admCoefficientManager::readCoefficient", dict)
            << "Unrecognized format looking up " << lookupName << " for "
            << coeffName << " coefficient."
            << abort(FatalIOError);
    }
    return * coeffPtr;
}


/*
void Foam::admCoefficientManager::saveState(const label slot)
{
    for
    (
        HashPtrTable<admCoefficient>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        iter()->saveState(slot);
    }
}


void Foam::admCoefficientManager::clearState(const label slot)
{
    for
    (
        HashPtrTable<admCoefficient>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        iter()->clearState(slot);
    }
}


void Foam::admCoefficientManager::loadState(const label slot)
{
    for
    (
        HashPtrTable<admCoefficient>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        iter()->loadState(slot);
    }
}


Foam::label Foam::admCoefficientManager::nStates() const
{
    HashPtrTable<admCoefficient>::const_iterator iter = begin();
    if (iter != end())
    {
        return iter()->nStates();
    }
    return 1;
}


bool Foam::admCoefficientManager::validState(const label slot) const
{
    HashPtrTable<admCoefficient>::const_iterator iter = begin();
    if (iter != end())
    {
        return iter()->validState(slot);
    }
    return true;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
