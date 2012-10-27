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

#include "admReactionManager.H"
//#include "admReaction.H"
#include "IOstreams.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admReactionManager, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admReactionManager::readDict()
{
//    rateInhibNames_ = admInhibitionDict_.toc();
    wordList rateInhibNames(admInhibitionDict_.toc());

    forAll(rateInhibNames, in)
    {
        autoPtr<admRateInhibition> tAdmRateInhib
        (
            admRateInhibition::New
            (
                runTime_,
                admVars_,
                admCoeffs_,
                admInhibitionDict_,
                rateInhibNames[in]
            )
        );

        inhibitions_.insert
        (
            rateInhibNames[in],
            tAdmRateInhib.ptr()
        );
    }

    wordList reacNames(admReactionDict_.toc());

    setSize(reacNames.size());
    forAll(reacNames, rn)
    {
        set
        (
            rn,
            new admReaction
            (
                admVars_,
                admCoeffs_,
                * this,
                admVars_.standard().size(),
                admVars_.implicit().size(),
                admVars_.derived().size(),
                admReactionDict_,
                reacNames[rn]
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admReactionManager::admReactionManager
(
    admTime& runTime,
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    const word& inhibitionDictName,
    const word& reactionDictName
)
:
    PtrList<admReaction>(0),
    runTime_(runTime),
    admVars_(admVars),
    admCoeffs_(admCoeffs),
    
    admInhibitionDict_
    (
        IOobject
        (
            inhibitionDictName,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    admReactionDict_
    (
        IOobject
        (
            reactionDictName,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    inhibitions_(100),
    
    initialized_(false)
{
    readDict();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::admReactionManager::initialize()
{
    if (initialized_)
    {
        FatalErrorIn("admReactionManager::initialize")
            << "Object already initialized."
            << abort(FatalError);
    }

    for (label reacIndex(0); reacIndex < size(); reacIndex++)
    {
        operator[](reacIndex).initialize();
    }
    
    initialized_ = true;
}


Foam::wordList Foam::admReactionManager::inhibitionToc() const
{
    return inhibitions_.toc();

/* When inhibitions_ was a PtrList
    wordList returnMe(inhibitions_.size());
    forAll(inhibitions_, inhibIndex)
    {
        returnMe[inhibIndex] = inhibitions_[inhibIndex].name();
    }
    return returnMe;*/
}


bool Foam::admReactionManager::inhibitionFound
(
    const word& inhibitionName
) const
{
    const admRateInhibition * inhibPtr(inhibitions_[inhibitionName]);
    return bool(inhibPtr);
    
/* When inhibitions_ was a PtrList
    forAll(inhibitions_, inhibIndex)
    {
        if (inhibitionName == inhibitions_[inhibIndex].name())
        {
            return true;
        }
    }
    return false;*/
}


const Foam::admRateInhibition& Foam::admReactionManager::inhibitionLookup
(
    const word& inhibitionName
) const
{
    const admRateInhibition * inhibPtr(inhibitions_[inhibitionName]);
#   ifdef FULLDEBUG
    if (!inhibPtr)
    {
        FatalErrorIn("admReactionManager::rateInhibition")
            << inhibitionName << " is not a valid inhibition name.  Valid "
            << "names are: " << inhibitions_.toc()
            << abort(FatalError);
    }
#   endif
    return * inhibPtr;
    
/* When inhibitions_ was a PtrList
    forAll(inhibitions_, inhibIndex)
    {
        if (inhibitionName == inhibitions_[inhibIndex].name())
        {
            return inhibitions_[inhibIndex];
        }
    }
    FatalErrorIn("admReactionManager::inhibitionLookup")
        << "admRateInhibition " << inhibitionName << " not found."
        << abort(FatalError);

    // dummy return
    return inhibitions_[0];*/
}


Foam::wordList Foam::admReactionManager::toc() const
{
    wordList returnMe(size());
    for (label reacIndex(0); reacIndex < size(); reacIndex++)
    {
        returnMe[reacIndex] = operator[](reacIndex).name();
    }
    return returnMe;
}


bool Foam::admReactionManager::found(const word& reactionName) const
{
    for (label reacIndex(0); reacIndex < size(); reacIndex++)
    {
        if (reactionName == operator[](reacIndex).name())
        {
            return true;
        }
    }
    return false;
}


const Foam::admReaction& Foam::admReactionManager::lookup
(
    const word& reactionName
) const
{
    for (label reacIndex(0); reacIndex < size(); reacIndex++)
    {
        if (reactionName == operator[](reacIndex).name())
        {
            return operator[](reacIndex);
        }
    }
    FatalErrorIn("admReactionManager::lookup")
        << "admReaction " << reactionName << " not found."
        << abort(FatalError);

    // dummy return
    return operator[](0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
