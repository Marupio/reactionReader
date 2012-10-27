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

//#include "admReactionManager.H"
#include "admReaction.H"
#include "admCoefficientManager.H"
//#include "admReactionRate.H"
#include "IOstreams.H"
#include "token.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admReaction, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admReaction::readDict(const dictionary& dict)
{
    yields_.setSize(ratePtr_().admVars().nAll());
    forAll(yields_, varIndex)
    {
        yields_.set
        (
            varIndex,
            &ratePtr_().admCoeffs().zero()
        );
    }

    const dictionary& yieldDict(dict.subDict(name_).subDict("yields"));
    wordList yieldVars(yieldDict.toc());
    forAll(yieldVars, i)
    {
        if (!ratePtr_().admVars().found(yieldVars[i]))
        {
            WarningIn("admReaction::readDict")
                << "Reaction " << name_ << " has yield defined for "
                << "unknown variable " << yieldVars[i] << ".  "
                << "Ignoring entry." << endl;
            continue;
        }

        // Index where the coefficient will go in yields_ list is:
        label varIndex
        (
            ratePtr_().admVars().lookup(yieldVars[i]).globalIndex()
        );
        
        yields_.set
        (
            varIndex,
            & ratePtr_().admCoeffs().readCoefficient
            (
                yieldDict,
                yieldVars[i],
                word(name_ + "(yield(" + yieldVars[i] + "))")
            )
        );
        
        /*
        // Pointer to coefficient that will go there:
        const admCoefficient * coeffPtr;
        
        // Check dictionary format - can be any of:
        //      varName     {subDict}  (full inline coefficient definition)
        //      varName     word; (gives coeffName)
        //      varName     dimensionedScalar; (short for constant coefficient)
        if (yieldDict.isDict(yieldVars[i]))
        {
            // Full inline coefficient definition
            // Create a new coefficient based on this dictionary
            word coeffName(name_ + "(yield(" + yieldVars[i] + "))");
            dictionary tempDict;
            tempDict.set(coeffName, yieldDict.subDict(yieldVars[i]));
            rate_.admCoeffs().addNew
            (
                tempDict,
                coeffName
            );
            coeffPtr = & rate_.admCoeffs()(coeffName);
        }
        else if
        (
            equationReader::isWord
            (
                yieldDict.lookup(yieldVars[i])
            )
        )
        {
            // Names a coefficient
            word coeffName(yieldDict.lookup(yieldVars[i]));
            coeffPtr = & rate_.admCoeffs()(coeffName);
        }
        else if
        (
            equationReader::isEquation
            (
                yieldDict.lookup(yieldVars[i])
            )
        )
        {
            // Shorthand for a new constant coefficient
            // Create a new constant coefficient based on entry
            word coeffName(name_ + "(yield(" + yieldVars[i] + "))");
            dictionary tempDict;
            tempDict.set(coeffName, yieldDict.lookup(yieldVars[i]));
            rate_.admCoeffs().addNew
            (
                tempDict,
                coeffName
            );
            coeffPtr = & rate_.admCoeffs()(coeffName);
        }
        
        yields_.set
        (
            varIndex,
            coeffPtr
        );*/
    } // end forAll yield vars
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admReaction::admReaction
(
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    admReactionManager& admReacs,
    label standardSize,
    label implicitSize,
    label derivedSize,
    const dictionary& dict,
    const word& reactionName
)
:
    runTime_(admVars.runTime()),
    mesh_(admVars.mesh()),

    name_(reactionName),
    
    ratePtr_
    (
        admReactionRate::New
        (
            runTime_,
            admVars,
            admCoeffs,
            admReacs,
            dict,
            reactionName
        )
    ),
//    rate_(ratePtr_()),
    yields_(0),
    nonZeroYields_(0),
    
    initialized_(false)
{
    readDict(dict);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::admReaction::initialize()
{
    if (initialized_)
    {
        FatalErrorIn("admReaction::initialize")
            << "Object already initialized."
            << abort(FatalError);
    }

    nonZeroYields_.clear();
    
    forAll(yields_, globalIndex)
    {
        if (yields_[globalIndex].evaluateNonZero())
        {
            label newIndex(nonZeroYields_.size());
            nonZeroYields_.setSize(newIndex + 1);
            nonZeroYields_[newIndex] = globalIndex;
        }
    }
    ratePtr_().initialize();
    initialized_ = true;
}

// ************************************************************************* //
