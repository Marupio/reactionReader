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

#include "admReactionRate.H"
#include "admReactionManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admReactionRate, 0);

    defineRunTimeSelectionTable(admReactionRate, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admReactionRate::readDict
(
    const dictionary& dict,
    const word& reactionName
)
{
    if (dict.subDict(reactionName).found("rateInhibitions"))
    {
        wordList rateInhibitionNames
        (
            dict.subDict(reactionName).lookup("rateInhibitions")
        );
        forAll(rateInhibitionNames, i)
        {
            if (admReacs_.inhibitionFound(rateInhibitionNames[i]))
            {
                label newIndex(rateInhibitions_.size());
                rateInhibitions_.setSize(newIndex + 1);
                rateInhibitions_.set
                (
                    newIndex,
                    &admReacs_.inhibitionLookup(rateInhibitionNames[i])
                );
            }
            else
            {
                WarningIn("admReactionRate::readDict")
                    << "Reaction " << reactionName << " lists unknown rate "
                    << "inhibition " << rateInhibitionNames[i] << ". Entry is "
                    << "ignored." << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admReactionRate::admReactionRate
(
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    admReactionManager& admReacs,
    const dictionary& dict,
    const word& name
)
:
    admCalculusInterface
    (
        admVars,
        "rate(" + name + ")",
        dimless // dimensions are expected to be set by derived classes
    ),

    admCoeffs_(admCoeffs),
    admReacs_(admReacs),

    rateInhibitions_(0),
    nonZeroInhibitionDdy_(0),
    initialized_(false),
        
    lastEvaluateValue_(mesh_.nCells()),
    lastInhibitionValue_(mesh_.nCells()),

    lastEvaluateMilestone_(-1),
    lastInhibitionMilestone_(-1)
{
    readDict(dict, name);
}


/*Foam::admReactionRate::admReactionRate(const admReactionRate& admc)
:
    admVars_(admc.admVars_),
    admCoeffs_(admc.admCoeffs_),
    admReacs_(admc.admReacs_),
    name_(admc.name_),
    dimensions_(admc.dimensions_),
    rateInhibitionNames_(admc.rateInhibitionNames_)
{}*/

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::admReactionRate::initialize()
{
    if (initialized_)
    {
        FatalErrorIn("admReactionRate::initialize")
            << "Object already initialized."
            << abort(FatalError);
    }

    nonZeroInhibitionDdy_.setSize(admVars_.nAll());

    forAll(admVars_.all(), globalIndex)
    {
        const admVariable& wrtVar(admVars_.all()[globalIndex]);
        forAll(rateInhibitions_, i)
        {
            bool nonZero(true);
            
            // d(I1 I2 I3 I4) / dy = ... + I1 d(I2)/dy I3 I4 + ...
            //                         Front^ ^Middle^ ^ Tail

            // Front
            for (label j(0); j < i; j++)
            {
                if (!rateInhibitions_[j].evaluateNonZero())
                {
                    nonZero = false;
                }
            }

            // Middle
            if (!rateInhibitions_[i].ddyNonZero(wrtVar))
            {
                nonZero = false;
            }

            // Tail
            for (label j(i + 1); j < rateInhibitions_.size(); j++)
            {
                if (!rateInhibitions_[j].evaluateNonZero())
                {
                    nonZero = false;
                }
            }
            if (nonZero)
            {
                label newIndex(nonZeroInhibitionDdy_[globalIndex].size());
                nonZeroInhibitionDdy_[globalIndex].setSize(newIndex + 1);
                nonZeroInhibitionDdy_[globalIndex][newIndex] = i;
            }
        }
    }
    initialized_ = true;
}


Foam::scalar Foam::admReactionRate::inhibition(label cellIndex) const
{
    scalar returnMe(1.0);
    if (lastInhibitionMilestone_ == runTime_.milestone())
    {
        return lastInhibitionValue_[cellIndex];
    }
    forAll(rateInhibitions_, inhibIndex)
    {
        returnMe *= rateInhibitions_[inhibIndex].evaluate(cellIndex);
    }
    return returnMe;
}


Foam::tmp<scalarField> Foam::admReactionRate::inhibitionField() const
{
    if (lastInhibitionMilestone_ != runTime_.milestone())
    {
        lastInhibitionValue_ = 1.0;
        forAll(rateInhibitions_, inhibIndex)
        {
            lastInhibitionValue_ *=
                rateInhibitions_[inhibIndex].evaluateField();
        }
        lastInhibitionMilestone_ = runTime_.milestone();
    }
    return tmp<scalarField>
    (
        new scalarField(lastInhibitionValue_)
    );
}


Foam::tmp<dimensionedScalarField>
    Foam::admReactionRate::inhibitionDimensionedField() const
{
    return tmp<dimensionedScalarField>
    (
        new dimensionedScalarField
        (
            IOobject
            (
                inhibitionName(),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            inhibitionDims(),
            inhibitionField()
        )
    );
}


bool Foam::admReactionRate::inhibitionNonZero() const
{
    forAll(rateInhibitions_, inhibIndex)
    {
        if (!rateInhibitions_[inhibIndex].evaluateNonZero())
        {
            return false;
        }
    }
    return true;
}


Foam::scalar Foam::admReactionRate::inhibitionDdy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // This is the chain rule applied to the inhibitions
    scalar returnMe(0.0);

    const label& wrtGlobalI(wrtVar.globalIndex());
    forAll(nonZeroInhibitionDdy_[wrtGlobalI], nzI)
    {
        label i(nonZeroInhibitionDdy_[wrtGlobalI][nzI]);
 //}The 4 lines above mimic:
 // forAll(rateInhibitions_, i)
 // {
    
        scalar product(1.0);
        
        // d(I1 I2 I3 I4) / dy = ... + I1 d(I2)/dy I3 I4 + ...
        //                         Front^ ^Middle^ ^ Tail

        // Front
        for (label j(0); j < i; j++)
        {
            product *= rateInhibitions_[j].evaluate(cellIndex);
        }

        // Middle
        product *= rateInhibitions_[i].ddy(wrtVar, cellIndex);

        // Tail
        for (label j(i + 1); j < rateInhibitions_.size(); j++)
        {
            product *= rateInhibitions_[j].evaluate(cellIndex);
        }
        
        returnMe += product;
    }
    return returnMe;
}


Foam::tmp<scalarField> Foam::admReactionRate::inhibitionDdyField
(
    const admVariable& wrtVar
) const
{
    // This is the chain rule applied to the inhibitions
    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            mesh_.nCells(),
            0.0
        )
    );
    scalarField& returnMe(tReturnMe());

    const label& wrtGlobalI(wrtVar.globalIndex());
    forAll(nonZeroInhibitionDdy_[wrtGlobalI], nzI)
    {
        label i(nonZeroInhibitionDdy_[wrtGlobalI][nzI]);
 //}The 4 lines above mimic:
 // forAll(rateInhibitions_, i)
 // {
    
        scalarField product(mesh_.nCells(), 1.0);
        
        // d(I1 I2 I3 I4) / dy = ... + I1 d(I2)/dy I3 I4 + ...
        //                         Front^ ^Middle^ ^ Tail

        // Front
        for (label j(0); j < i; j++)
        {
            product *= rateInhibitions_[j].evaluateField();
        }

        // Middle
        product *= rateInhibitions_[i].ddyField(wrtVar);

        // Tail
        for (label j(i + 1); j < rateInhibitions_.size(); j++)
        {
            product *= rateInhibitions_[j].evaluateField();
        }
        
        returnMe += product;
    }
    return tReturnMe;
}


Foam::tmp<dimensionedScalarField>
    Foam::admReactionRate::inhibitionDdyDimensionedField
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
                inhibitionDdyName(wrtVar),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            inhibitionDdyDims(wrtVar),
            inhibitionDdyField(wrtVar)
        )
    );
}


bool Foam::admReactionRate::inhibitionDdyNonZero
(
    const admVariable& wrtVar
) const
{
    return (nonZeroInhibitionDdy_[wrtVar.globalIndex()].size() != 0);
}


Foam::tmp<dimensionedScalarField>
    Foam::admReactionRate::uninhibitedDimensionedField() const
{
    return tmp<dimensionedScalarField>
    (
        new dimensionedScalarField
        (
            IOobject
            (
                uninhibitedName(),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            uninhibitedDims(),
            uninhibitedField()
        )
    );
}


Foam::tmp<dimensionedScalarField>
    Foam::admReactionRate::uninhibitedDdyDimensionedField
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
                uninhibitedDdyName(wrtVar),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            uninhibitedDdyDims(wrtVar),
            uninhibitedDdyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admReactionRate::evaluate(label cellIndex) const
{
    if (lastEvaluateMilestone_ != runTime_.milestone())
    {
        return uninhibited(cellIndex) * inhibition(cellIndex);
    }
    return lastEvaluateValue_[cellIndex];
}


Foam::tmp<scalarField> Foam::admReactionRate::evaluateField() const
{

    if (lastEvaluateMilestone_ != runTime_.milestone())
    {
        lastEvaluateValue_ = uninhibitedField() * inhibitionField();
        lastEvaluateMilestone_ = runTime_.milestone();
    }
    return tmp<scalarField>
    (
        new scalarField(lastEvaluateValue_)
    );
}


bool Foam::admReactionRate::evaluateNonZero() const
{
    return true;
}


Foam::scalar Foam::admReactionRate::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return uninhibited(cellIndex) * inhibitionDdy(wrtVar, cellIndex) +
        uninhibitedDdy(wrtVar, cellIndex) * inhibition(cellIndex);
}


Foam::tmp<scalarField> Foam::admReactionRate::ddyField
(
    const admVariable& wrtVar
) const
{
    return uninhibitedField() * inhibitionDdyField(wrtVar) +
        uninhibitedDdyField(wrtVar) * inhibitionField();
}


bool Foam::admReactionRate::ddyNonZero(const admVariable& wrtVar) const
{
    return (uninhibitedDdyNonZero(wrtVar) || inhibitionDdyNonZero(wrtVar));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newAdmReactionRate.C"

// ************************************************************************* //
