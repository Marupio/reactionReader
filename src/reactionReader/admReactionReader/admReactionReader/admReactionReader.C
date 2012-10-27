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

#include "admReactionReader.H"
#include "IOstreams.H"
#include "token.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admReactionReader, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admReactionReader::initialize()
{
    if (initialized_)
    {
        FatalErrorIn("admReactionReader::initialize")
            << "Object already initialized.  This is considered fatal."
            << abort(FatalError);
    }
    //admVars_.initialize();    // triggered by admCoeffs constructor
    //admCoeffs_.initialize();  // no initialization routine for coeffs
    admReacs_.initialize();
    createNonZeroYields();
    initialized_ = true;
}


void Foam::admReactionReader::createNonZeroYields()
{
    nonZeroYields_.clear();
    nonZeroStandardYields_.clear();
    nonZeroImplicitYields_.clear();
    nonZeroDerivedYields_.clear();
    nonZeroYields_.setSize(admVars_.nAll());
    nonZeroStandardYields_.setSize(admVars_.nStandard());
    nonZeroImplicitYields_.setSize(admVars_.nImplicit());
    nonZeroDerivedYields_.setSize(admVars_.nDerived());

    // All variables
    forAll(admVars_.all(), globalIndex)
    {
        const admVariable& theVariable(admVars_.all(globalIndex));
        nonZeroYields_.set
        (
            globalIndex,
            new UPtrList<const admReaction>
        );
        forAll(admReacs_, reacIndex)
        {
            const admReaction& theReaction(admReacs_[reacIndex]);
            if (theReaction.yield(theVariable).evaluateNonZero())
            {
                label newIndex(nonZeroYields_[globalIndex].size());
                nonZeroYields_[globalIndex].setSize(newIndex + 1);
                nonZeroYields_[globalIndex].set
                (
                    newIndex,
                    & theReaction
                );
            }
        }
    }

    // Standard variables
    forAll(admVars_.standard(), localIndex)
    {
        const admStandardVariable& theVariable(admVars_.standard(localIndex));
        nonZeroStandardYields_.set
        (
            localIndex,
            new UPtrList<const admReaction>
        );
        forAll(admReacs_, reacIndex)
        {
            const admReaction& theReaction(admReacs_[reacIndex]);
            if (theReaction.yield(theVariable).evaluateNonZero())
            {
                label newIndex(nonZeroStandardYields_[localIndex].size());
                nonZeroStandardYields_[localIndex].setSize(newIndex + 1);
                nonZeroStandardYields_[localIndex].set
                (
                    newIndex,
                    & theReaction
                );
            }
        }
    }

    // Implicit variables
    forAll(admVars_.implicit(), localIndex)
    {
        const admImplicitVariable& theVariable(admVars_.implicit(localIndex));
        nonZeroImplicitYields_.set
        (
            localIndex,
            new UPtrList<const admReaction>
        );
        forAll(admReacs_, reacIndex)
        {
            const admReaction& theReaction(admReacs_[reacIndex]);
            if (theReaction.yield(theVariable).evaluateNonZero())
            {
                label newIndex(nonZeroImplicitYields_[localIndex].size());
                nonZeroImplicitYields_[localIndex].setSize(newIndex + 1);
                nonZeroImplicitYields_[localIndex].set
                (
                    newIndex,
                    & theReaction
                );
            }
        }
    }

    // Derived variables
    forAll(admVars_.derived(), localIndex)
    {
        const admDerivedVariable& theVariable(admVars_.derived(localIndex));
        nonZeroDerivedYields_.set
        (
            localIndex,
            new UPtrList<const admReaction>
        );
        forAll(admReacs_, reacIndex)
        {
            const admReaction& theReaction(admReacs_[reacIndex]);
            if (theReaction.yield(theVariable).evaluateNonZero())
            {
                label newIndex(nonZeroDerivedYields_[localIndex].size());
                nonZeroDerivedYields_[localIndex].setSize(newIndex + 1);
                nonZeroDerivedYields_[localIndex].set
                (
                    newIndex,
                    & theReaction
                );
            }
        } // end forAll reactions
    } // end forAll derived
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admReactionReader::admReactionReader
(
    admTime& runTime,
    fvMesh& mesh,
    const word& variableDictName,
    const word& coefficientDictName,
    const word& inhibitionDictName,
    const word& reactionDictName,
    const word& settingsDictName,
    const word& hooksDictName
)
:
    regIOobject
    (
        IOobject
        (
            "admReactionReader",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    
    initialized_(false),
    
    runTime_(runTime),
    mesh_(mesh),

    admSettingsDict_
    (
        IOobject
        (
            settingsDictName,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    admEquationDict_
    (
        runTime_.eqnDict()
    ),

    admVars_
    (
        *this,
        variableDictName
    ),

    admCoeffs_
    (
        runTime_,
        mesh_,
        admVars_,
        coefficientDictName
    ),

    admReacs_
    (
        runTime_,
        admVars_,
        admCoeffs_,
        inhibitionDictName,
        reactionDictName
    ),
    
    outputFlagsDict_
    (
        admSettingsDict_.found("outputFlags")
      ? admSettingsDict_.subDict("outputFlags")
      : dictionary::null
    ),
    outputTimeDetails_
    (
        outputFlagsDict_.found("timeDetails")
      ? bool(Switch(outputFlagsDict_.lookup("timeDetails")))
      : false
    ),
    outputReactionAverages_
    (
        outputFlagsDict_.found("reactionAverages")
      ? bool(Switch(outputFlagsDict_.lookup("reactionAverages")))
      : true
    )
{
    // Initialize model components
    initialize();

    // Override admTime::debug
    if (outputTimeDetails_)
    {
        admTime::debug = 1;
    }
    else
    {
        admTime::debug = 0;
    }
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::admReactionReader::~admReactionReader()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admReactionReader::calculateRms
(
    const scalarField& sfIn
)
{
    if (sfIn.size())
    {
//        scalar rms(0.0);
//        TFOR_ALL_S_OP_FUNC_F(scalar, rms, +=, sqr, scalar, sfIn)
//        rms /= sfIn.size();
//        return sqrt(rms);
        
        return sqrt(gSumSqr(sfIn) / sfIn.size());
    }
    else
    {
        WarningIn("admReactionReader::calculateRms")
            << "Field is empty, assuming rms = 0"
            << endl;
        return scalar(0.0);
    }
}


Foam::scalar Foam::admReactionReader::calculateRms
(
    const scalarField& sfIn,
    const scalarField& weights
)
{
    scalarField weighted(sfIn * weights);

    if (sfIn.size())
    {
        return sqrt(gSumSqr(weighted) / gSum(weights));
    }
    else
    {
        WarningIn("admReactionReader::calculateRms")
            << "Field is empty, assuming rms = 0"
            << endl;
        return scalar(0.0);
    }
}


Foam::scalar Foam::admReactionReader::calculateRmsError
(
    const scalarField& a,
    const scalarField& b,
    scalar scalingFactor
)
{
    return calculateRms(a - b) / scalingFactor;
}


Foam::scalar Foam::admReactionReader::calculateRmsError
(
    const scalarField& a,
    const scalarField& b,
    const scalarField& weights,
    scalar scalingFactor
)
{
    return calculateRms(a - b, weights) / scalingFactor;
}


Foam::scalar Foam::admReactionReader::calculateNextDeltaT
(
    const scalarField& coarseResult,
    const scalarField& fineResult,
    const scalar& convergence,
    const scalarField& weights,
    scalar currentDeltaT,
    scalar scalingFactor
)
{
    scalar rmsError
    (
        stabilise
        (
            calculateRmsError
            (
                coarseResult,
                fineResult,
                weights,
                scalingFactor
            ),
            VSMALL
        )
    );
    
    return currentDeltaT * pow(convergence / rmsError, 0.2);
}


void Foam::admReactionReader::setMilestone()
{
    runTime_.setMilestone();
    admVars_.correctAllBoundaryConditions();
}


void Foam::admReactionReader::saveState(const label slot)
{
    if ((slot > nStates()) || (slot < 0))
    {
        FatalErrorIn("admReactionReader::saveState")
            << "saveState slot " << slot << " out of range: 0, "
            << nStates() << ". Can only save over existing "
            << "slots or create one at end."
            << abort(FatalError);
    }
    if (slot == nStates())
    {
        validStates_.setSize(slot + 1);
    }
    validStates_[slot] = true;
    runTime_.saveState(slot);
    admVars_.saveState(slot);
    //admCoeffs_.saveState(slot);
    //admReacs_.saveState(slot);
}


void Foam::admReactionReader::clearState(const label slot)
{
    if ((slot >= nStates()) || (slot < 0))
    {
        WarningIn("admReactionReader::clearState")
            << "Attempting to clear non-existent slot " << slot << ". "
            << "Ignoring. Range = 0, " << nStates() << endl;
        return;
    }
    validStates_[slot] = false;
    runTime_.clearState(slot);
    admVars_.clearState(slot);
    //admCoeffs_.clearState(slot);
    //admReacs_.clearState(slot);
}


void Foam::admReactionReader::loadState(const label slot)
{
    if ((slot >= nStates()) || (slot < 0))
    {
        FatalErrorIn("admReactionReader::loadState")
            << "Load slot " << slot << " out of range 0, " << nStates()
            << abort(FatalError);
    }
    if (!validStates_[slot])
    {
        FatalErrorIn("admModhel::loadState")
            << "Attempting to load empty save state " << slot
            << abort(FatalError);
    }
    runTime_.loadState(slot);
    admVars_.loadState(slot);
    //admCoeffs_.loadState(slot);
    //admReacs_.loadState(slot);
}


Foam::label Foam::admReactionReader::nStates() const
{
    return validStates_.size();
}


bool Foam::admReactionReader::validState(const label slot) const
{
    // Is slot out-of-range?
    if ((slot >= nStates()) || (slot < 0))
    {
        return false;
    }

    return validStates_[slot];
}


bool Foam::admReactionReader::writeData(Ostream&) const
{
    // Do nothing
    return true;
}


void Foam::admReactionReader::reportVariablesHeader(Ostream& os) const
{
    os << "Time,Cell";
    forAll(admVars_.all(), varIndex)
    {
        os << "," << admVars_.all(varIndex).name();
    }
    os << endl;
}


void Foam::admReactionReader::reportVariables(Ostream& os, const label cellIndex) const
{
    os << runTime_.value() << "," << cellIndex;
    forAll(admVars_.all(), varIndex)
    {
        os << "," << admVars_.all(varIndex).evaluate(cellIndex);
    }
    os << endl;
}


void Foam::admReactionReader::reportReactionRatesHeader(Ostream& os) const
{
    os << "Time,Cell";
    forAll(admReacs_, reacI)
    {
        os << "," << admReacs_[reacI].name();
    }
    os << endl;
}


void Foam::admReactionReader::reportReactionRates
(
    Ostream& os,
    const label cellIndex
) const
{
    os << runTime_.value() << "," << cellIndex;
    forAll(admReacs_, reacI)
    {
        os << "," << admReacs_[reacI].rate().evaluate(cellIndex);
    }
    os << endl;
}


void Foam::admReactionReader::reportCoefficientsHeader(Ostream& os) const
{
    os << "Time,Cell";
    for
    (
        HashPtrTable<admCoefficient>::const_iterator iter = admCoeffs_.begin();
        iter != admCoeffs_.end();
        ++iter
    )
    {
        os << "," << iter()->name();
    }
    os << endl;
}


void Foam::admReactionReader::reportCoefficients
(
    Ostream& os,
    const label cellIndex
) const
{
    os << runTime_.value() << "," << cellIndex;
    for
    (
        HashPtrTable<admCoefficient>::const_iterator iter = admCoeffs_.begin();
        iter != admCoeffs_.end();
        ++iter
    )
    {
        os << "," << iter()->evaluate(cellIndex);
    }
    os << endl;
}

// ************************************************************************* //
