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

#include "admTime.H"
#include "equation.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admTime, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from objectRegistry arguments
Foam::admTime::admTime
(
    const word& name,
    const fileName& rootPath,
    const fileName& caseName,
    const fileName& systemName,
    const fileName& constantName,
    const word& equationsDictName
)
:
    Time
    (
        name,
        rootPath,
        caseName,
        systemName,
        constantName
    ),
    
    eqns_
    (
        IOobject
        (
            "equations",
            timeName(),
            * this,
            IOobject::NO_READ,
            controlDict().found("outputEquations")
          ? Switch(controlDict().lookup("outputEquations")) == true
            ? IOobject::AUTO_WRITE
            : IOobject::NO_WRITE
          : IOobject::NO_WRITE
        ),        
        controlDict().found("outputEquationDataSources")
      ? bool(Switch(controlDict().lookup("outputEquationDataSources")))
      : false
    ),
    
    admEquationDict_
    (
        IOobject
        (
            equationsDictName,
            constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    doubled_(false),
    
    maxDeltaT_
    (
        controlDict().found("maxDeltaT")
      ? readScalar(controlDict().lookup("maxDeltaT"))
      : VGREAT
    ),
    
    minDeltaT_
    (
        controlDict().found("minDeltaT")
      ? readScalar(controlDict().lookup("minDeltaT"))
      : 0.0
    ),
    
    lastOutputTime_(startTime_),
    minimumOutputSpacing_
    (
        controlDict().found("minimumOutputSpacing")
      ? readScalar(controlDict().lookup("minimumOutputSpacing"))
      : 0.0
    ),
    
    timeValue_a(0),
    timeIndex_a(0),
    deltaT_a(0),
    deltaTSave_a(0),
    deltaT0_a(0),
    deltaTchanged_a(0),
    outputTimeIndex_a(0),
    outputTime_a(0),
    doubled_a(0)
{
    if (controlDict().found("resumeDeltaT"))
    {
        Switch resumeDeltaT(controlDict().lookup("resumeDeltaT"));
        if (resumeDeltaT && (deltaT0_ > 0))
        {
            deltaT_ = deltaT0_;
        }
    }
    
    eqns_.addSource(admEquationDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::admTime::readScalarOrEquation
(
    Istream& is,
    const word& varName
)
{
    dimensionedScalar ds(readDimensionedScalarOrEquation(is, varName));
    if (!ds.dimensions().dimensionless())
    {
        WarningIn("admTime::readScalarOrEquation")
            << "Reading " << varName << " from dictionary.  Expecting "
            << "dimensionless, received: " << ds.dimensions() << ". Ignoring "
            << "dimensions." << endl;
    }
    return ds.value();
}


Foam::dimensionedScalar Foam::admTime::readDimensionedScalarOrEquation
(
    Istream& is,
    const word& varName
)
{
    dimensionedScalar ds
    (
        varName,
        dimless,
        0.0
    );
    bool hasDimensions(false);

    tokenList tl(13);
    label found(0);
    while (!is.eof())
    {
        tl[found] = token(is);
        found++;
        if (found > 12)
        {
            FatalIOErrorIn("admTime::readDimensionedScalarOrEquation", is)
                << varName << " has too many tokens (bad format).  Valid "
                << "formats are:" << token::NL << token::TAB
                << "keyword     scalar;" << token::NL << token::TAB
                << "keyword     \"equation\";" << token::NL << token::TAB
                << "keyword     [dimensionSet] scalar;"
                << token::NL << token::TAB
                << "keyword     [dimensionSet] \"equation\";"
                << token::NL << token::TAB
                << "keyword     ignoredWord [dimensionSet] scalar;"
                << token::NL << token::TAB
                << "keyword     ignoredWord [dimensionSet] \"equation\";"
                << exit(FatalIOError);
        }
    }
    tl.setSize(found);
    is.rewind();
    label i(0);
    
    if (tl[i].isWord())
    {
        // Starts with optional ignoredWord
        word ignoredWord(is);
        i++;
    }
    if (tl[i].isPunctuation())
    {
        // Has optional dimensionSet
        hasDimensions = true;
        ds.dimensions().reset(dimensionSet(is));
        i++;
        
        // End of dimensionSet is at next punctuation
        while (!tl[i].isPunctuation())
        {
            i++;
        }
        i++;
    }
    if (tl[i].isNumber())
    {
        ds.value() = readScalar(is);
        i++;
    }
    else if (tl[i].isString())
    {
        word eqnName("_once_" + varName);
        string rawText(tl[i].stringToken());
        equation eqn
        (
            eqnName,
            rawText,
            ds.dimensions(),
            hasDimensions
        );
        label eqnIndex(eqns_.readEquation(eqn));
        dimensionedScalar ds2
        (
            eqns_.evaluateDimensionedScalar(eqnIndex)
        );
        ds.value() = ds2.value();
        ds.dimensions().reset(ds2.dimensions());
        eqns_.deleteEquation(eqnIndex);
    }
    i++;
    if (i < found)
    {
        OStringStream badBit;
        for (; i < found; i++)
        {
            badBit << tl[i];
        }
        FatalIOErrorIn("admTime::readDimensionedScalarOrEquation", is)
            << "Invalid format reading " << varName << " from dictionary. "
            << "Extra tokens do not make sense:" << token::NL << token::TAB
            << badBit.str() << token::NL << "Valid formats are:" << token::NL
            << token::TAB << "keyword     scalar;" << token::NL << token::TAB
            << "keyword     \"equation\";" << token::NL << token::TAB
            << "keyword     [dimensionSet] scalar;"
            << token::NL << token::TAB
            << "keyword     [dimensionSet] \"equation\";"
            << token::NL << token::TAB
            << "keyword     ignoredWord [dimensionSet] scalar;"
            << token::NL << token::TAB
            << "keyword     ignoredWord [dimensionSet] \"equation\";"
            << exit(FatalIOError);
    }
    return ds;
}


void Foam::admTime::doubleDeltaT()
{
    deltaT_ += deltaT_;
    deltaTchanged_ = true;

    if (doubled_)
    {
        WarningIn("admTime::doubleDeltaT")
        << "Timestep has already been doubled." << endl;
    }
    doubled_ = true;
    
    if (debug)
    {
        Info << "admTime::doubleDeltaT : Doubling delta T to " << deltaT_
            << endl;
    }
}


void Foam::admTime::halveDeltaT()
{
    deltaT_ /= 2;
    deltaTchanged_ = true;
    
    if (!doubled_)
    {
        WarningIn("admTime::halveDeltaT")
        << "Timestep has already been halved." << endl;
    }
    doubled_ = false;
    
    if (debug)
    {
        Info << "admTime::halveDeltaT : Halving delta T to " << deltaT_
            << endl;
    }
}


void Foam::admTime::setDeltaTCoarse(const scalar deltaT)
{
    setDeltaTLimited(deltaT);
    if (writeControl_ == wcAdjustableRunTime)
    {
        scalar timeToNextWrite = max
        (
            0.0,
            (outputTimeIndex_ + 1) * writeInterval_
          - (value() - startTime_)
        );

        label nStepsToNextWrite
        (
            (timeToNextWrite / deltaT_ / 2 - SMALL) + 1
        );

        if (nStepsToNextWrite < 3)
        {
            deltaT_ = timeToNextWrite / nStepsToNextWrite / 2;
        }
    }
    
    if (debug)
    {
        Info << "admTime::setDeltaTCoarse : Setting delta T to " << deltaT_
            << endl;
    }
}


void Foam::admTime::setDeltaTLimited(const scalar deltaT)
{
    scalar limitedDeltaT(min(deltaT, maxDeltaT_));
    limitedDeltaT = max(limitedDeltaT, minDeltaT_);
    limitedDeltaT = min(limitedDeltaT, (endTime_ - this->value()) / 2);
    if (limitedDeltaT > 0)
    {
        deltaT_ = limitedDeltaT;
    }
    else
    {
        deltaT_ = deltaT;
    }
    deltaTchanged_ = true;

    if (debug)
    {
        Info << "admTime::setDeltaT : requested " << deltaT
            << ", limited to " << deltaT_ << endl;
    }
}


void Foam::admTime::setDeltaT(const scalar deltaT)
{
    deltaT_ = deltaT;
    deltaTchanged_ = true;
    if (debug)
    {
        Info << "admTime::setDeltaT : Setting delta T to " << deltaT_ << endl;
    }
}




void Foam::admTime::setTime
(
    const scalar& newTimeValue,
    const label& newTimeIndex,
    const scalar& newDeltaT,
    const scalar& newDeltaTSave,
    const scalar& newDeltaT0,
    const bool& newDeltaTchanged,
    const label& newOutputTimeIndex,
    const bool& newOutputTime
)
{
    Time::setTime(newTimeValue, newTimeIndex);
    deltaT_ = newDeltaT;
    deltaTSave_ = newDeltaTSave;
    deltaT0_ = newDeltaT0;
    deltaTchanged_ = newDeltaTchanged;
    outputTimeIndex_ = newOutputTimeIndex;
    outputTime_ = newOutputTime;

    if (debug)
    {
        Info << "admTime::setTime : Setting runTime to " << value()
            << " with dt = " << deltaT_ << endl;
    }
}


bool Foam::admTime::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    if (outputTime_)
    {
        lastOutputTime_ = value();
    }
    return Time::writeObject(fmt, ver, cmp);
}


void Foam::admTime::saveState(const label slot)
{
    if (debug)
    {
        Info << "admTime::saveState : Saving state " << slot << " with t = "
            << value() << ", dt = " << deltaT_ << endl;
    }
    if ((slot > nStates()) || (slot < 0))
    {
        FatalErrorIn("admTime::saveState")
            << "saveState slot " << slot << " out of range: 0, "
            << nStates() << ". Can only save over existing slots "
            << "or create one at end."
            << abort(FatalError);
    }
    if (slot == nStates())
    {
        timeValue_a.setSize(slot + 1);
        timeIndex_a.setSize(slot + 1);
        deltaT_a.setSize(slot + 1);
        deltaTSave_a.setSize(slot + 1);
        deltaT0_a.setSize(slot + 1);
        deltaTchanged_a.setSize(slot + 1);
        outputTimeIndex_a.setSize(slot + 1);
        outputTime_a.setSize(slot + 1);
        doubled_a.setSize(slot + 1);
    }
    timeValue_a[slot] = value();
    timeIndex_a[slot] = timeIndex_;
    deltaT_a[slot] = deltaT_;
    deltaTSave_a[slot] = deltaTSave_;
    deltaT0_a[slot] = deltaT0_;
    deltaTchanged_a[slot] = deltaTchanged_;
    outputTimeIndex_a[slot] = outputTimeIndex_;
    outputTime_a[slot] = outputTime_;
    doubled_a[slot] = doubled_;
}


void Foam::admTime::clearState(const label slot)
{
    if ((slot > nStates()) || (slot < 0))
    {
        WarningIn("admTime::clearState")
            << "Attempting to clear non-existent slot " << slot << ". "
            << "Ignoring. Range = 0, " << nStates() - 1 << endl;
        return;
    }
    timeValue_a[slot] = 0.0;
    timeIndex_a[slot] = -1;
    deltaT_a[slot] = 0.0;
    deltaTSave_a[slot] = 0.0;
    deltaT0_a[slot] = 0.0;
    deltaTchanged_a[slot] = false;
    outputTimeIndex_a[slot] = -1;
    outputTime_a[slot] = false;
    doubled_a[slot] = false;
}


void Foam::admTime::loadState(const label slot)
{
    if ((slot >= nStates()) || (slot < 0))
    {
        FatalErrorIn("admTime::loadState")
            << "Load slot " << slot << " out of range 0, " << nStates() - 1
            << abort(FatalError);
    }
    if (!validState(slot))
    {
        FatalErrorIn("admTime::loadState")
            << "Attempting to load empty save state " << slot
            << abort(FatalError);
    }
    setTime
    (
        timeValue_a[slot],
        timeIndex_a[slot],
        deltaT_a[slot],
        deltaTSave_a[slot],
        deltaT0_a[slot],
        deltaTchanged_a[slot],
        outputTimeIndex_a[slot],
        outputTime_a[slot]
    );
    doubled_ = doubled_a[slot];

    // Increment milestone - force all stored quantities to be recalculated
    setMilestone();

    if (debug)
    {
        Info << "admTime::loadState : Loading state " << slot << " with t = "
            << value() << ", dt = " << deltaT_ << endl;
    }
}


Foam::label Foam::admTime::nStates() const
{
    return timeValue_a.size();
}


bool Foam::admTime::validState(const label slot) const
{
    if ((slot >= timeValue_a.size()) || (slot < 0))
    {
        return false;
    }
    if (timeIndex_a[slot] == -1)
    {
        return false;
    }
    return true;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::admTime& Foam::admTime::plusPlusNoOutput()
{
    label oldOutputTimeIndex(outputTimeIndex_);
    Time::operator++();
    outputTimeIndex_ = oldOutputTimeIndex;

    if (debug)
    {
        Info << "admTime::plusPlusNoOutput : Incrementing runTime to "
            << value() << endl;
    }
    return * this;
}


Foam::admTime& Foam::admTime::operator++()
{
    // Check for timestep underrun
    if
    (
        (value() + deltaT().value()) == value()
    )
    {
        FatalErrorIn("admTime::operator++")
            << "Timestep underrun at " << timeName()
            << " + " << deltaT().value()
            << abort(FatalError);
    }

    Time::operator++();
    if
    (
        outputTime_
     && (
            (value() - lastOutputTime_) < minimumOutputSpacing_
        )
    )
    {
        outputTime_ = false;
    }
    if (!run() && !outputTime_)
    {
        outputTime_ = true;
    }
    
    if (debug)
    {
        Info << "admTime::operator++ : Incrementing runTime to " << value()
            << endl;
    }
    return * this;
}


Foam::admTime& Foam::admTime::operator++(int)
{
    return operator++();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
