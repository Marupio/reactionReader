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

#include "admVariableManager.H"
#include "IOstreams.H"
#include "token.H"
#include "admReactionReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admVariableManager, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admVariableManager::readDict()
{
    // Read defaults, if present
    if (admVariableDict_.found("defaults"))
    {
        dictionary& defaultsDict(admVariableDict_.subDict("defaults"));
        if (defaultsDict.found("standard"))
        {
            dictionary& standardDefaultsDict(defaultsDict.subDict("standard"));
            if (standardDefaultsDict.found("convergence"))
            {
                defaultStandardConvergence_ =
                    readScalar(standardDefaultsDict.lookup("convergence"));
                if (defaultStandardConvergence_ < 0)
                {
                    FatalIOErrorIn
                    (
                        "admVariableManager::readDict",
                        standardDefaultsDict
                    ) << "Default standard convergence cannot be negative"
                        << exit(FatalIOError);
                }
            }
            if (standardDefaultsDict.found("minErrorScale"))
            {
                defaultStandardMinErrorScale_ =
                    readScalar(standardDefaultsDict.lookup("minErrorScale"));
                if (defaultStandardMinErrorScale_ < 0)
                {
                    FatalIOErrorIn
                    (
                        "admVariableManager::readDict",
                        standardDefaultsDict
                    ) << "Default standard minErrorScale cannot be negative"
                        << exit(FatalIOError);
                }
            }
        }
        if (defaultsDict.found("implicit"))
        {
            dictionary& implicitDefaultsDict(defaultsDict.subDict("implicit"));
            if (implicitDefaultsDict.found("convergence"))
            {
                defaultImplicitConvergence_ =
                    readScalar(implicitDefaultsDict.lookup("convergence"));
                if (defaultImplicitConvergence_ < 0)
                {
                    FatalIOErrorIn
                    (
                        "admVariableManager::readDict",
                        implicitDefaultsDict
                    ) << "Default implicit convergence cannot be negative"
                        << exit(FatalIOError);
                }
            }
            if (implicitDefaultsDict.found("minErrorScale"))
            {
                defaultImplicitMinErrorScale_ =
                    readScalar(implicitDefaultsDict.lookup("minErrorScale"));
                if (defaultImplicitMinErrorScale_ < 0)
                {
                    FatalIOErrorIn
                    (
                        "admVariableManager::readDict",
                        implicitDefaultsDict
                    ) << "Default implicit minErrorScale cannot be negative"
                        << exit(FatalIOError);
                }
            }
        }
    }

    wordList varNames(admVariableDict_.toc());

    label autoSolveIndex(0);
    
    forAll(varNames, i)
    {
        if (varNames[i] == "defaults") continue;
        dictionary& varDict(admVariableDict_.subDict(varNames[i]));

        word varType;
        varType = word(varDict.lookup("type"));

        if (varType == "standard")
        {
            label localIndex(standard_.size());
            label globalIndex(all_.size());
            standard_.setSize(localIndex + 1);
            standard_.set
            (
                localIndex,
                new admStandardVariable
                (
                    runTime_,
                    * this,
                    varDict,
                    globalIndex,
                    localIndex,
                    varNames[i]
                )
            );
            
            all_.setSize(globalIndex + 1);
            all_.set
            (
                globalIndex,
                & standard_[localIndex]
            );
            
            // Add to changedByUdf or unchangedByUdf sublist
            if (standard_[localIndex].changedByUdf())
            {
                label newIndex(changedByUdf_.size());
                changedByUdf_.setSize(newIndex + 1);
                changedByUdf_.set(newIndex, &standard_[localIndex]);
            }
            else
            {
                label newIndex(unchangedByUdf_.size());
                unchangedByUdf_.setSize(newIndex + 1);
                unchangedByUdf_.set(newIndex, &standard_[localIndex]);
            }
        }
        else if (varType == "implicit")
        {
            label localIndex(implicit_.size());
            label globalIndex(all_.size());
            implicit_.setSize(localIndex + 1);
            implicit_.set
            (
                localIndex,
                new admImplicitVariable
                (
                    runTime_,
                    * this,
                    varDict,
                    globalIndex,
                    localIndex,
                    varNames[i]
                )
            );
            
            all_.setSize(globalIndex + 1);
            all_.set
            (
                globalIndex,
                & implicit_[localIndex]
            );
            
            // Add to implicitAutoSolve sublist if applicable
            if (implicit_[localIndex].autoSolve())
            {
                label newIndex(implicitAutoSolve_.size());
                implicitAutoSolve_.setSize(newIndex + 1);
                implicitAutoSolve_.set(newIndex, &implicit_[localIndex]);
                implicit_[localIndex].autoSolveIndex() = autoSolveIndex;
                autoSolveIndex++;
            }
        }
        else if (varType == "derived")
        {
            label localIndex(derived_.size());
            label globalIndex(all_.size());
            derived_.setSize(localIndex + 1);
            derived_.set
            (
                localIndex,
                new admDerivedVariable
                (
                    runTime_,
                    * this,
                    varDict,
                    globalIndex,
                    localIndex,
                    varNames[i]
                )
            );
            
            all_.setSize(globalIndex + 1);
            all_.set
            (
                globalIndex,
                & derived_[localIndex]
            );
        }
        else
        {
            // Fail
            FatalErrorIn("admVariableManager::readDict")
                << "Variable " << varNames[i] << " has unknown type."
                << " Available types are: standard, implicit, derived"
                << abort(FatalError);
        }
    } // forAll
    
    standard_.setSize(standard_.size());
    implicit_.setSize(implicit_.size());
    derived_.setSize(derived_.size());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admVariableManager::admVariableManager
(
    admReactionReader& model,
    const word& variableDictName
)
:
    model_(model),
    runTime_(model.runTime()),
    mesh_(model.mesh()),

    admVariableDict_
    (
        IOobject
        (
            variableDictName,
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    chainRuleSearchDepth_
    (
        readLabel(model.admSettingsDict().lookup("chainRuleSearchDepth"))
    ),

    defaultStandardConvergence_(-VGREAT),
    defaultImplicitConvergence_(-VGREAT),
    defaultStandardMinErrorScale_(-VGREAT),
    defaultImplicitMinErrorScale_(-VGREAT),

    nullStandard_(runTime_, *this),
    nullImplicit_(runTime_, *this),
    nullDerived_(runTime_, *this),

    standard_(0),
    implicit_(0),
    derived_(0),
    all_(0),
    
    changedByUdf_(0),
    unchangedByUdf_(0),
    implicitAutoSolve_(0),

    averagesDict_
    (
        IOobject
        (
            "averages",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),

    initialized_(false),

    standardRestore_(0),
    implicitRestore_(0)
{
    readDict();

    // Create Jacobian and reverse function indices (derived only)
    forAll(derived_, i)
    {
        derived_[i].createIndices();
    }

    // Create chainRule indices (derived only)
    if (chainRuleSearchDepth_ > 1)
    {
        forAll(derived_, i)
        {
            derived_[i].createChainRuleIndices(chainRuleSearchDepth_);
        }
    }

    // Add miscellaneous items to the equationReader
    //- mesh coordinates
    runTime_.eqns().vectorSources().addSource(mesh_.C());

    //- mesh volumes
    runTime_.eqns().scalarSources().addSource(mesh_.V());
    
    //- current time value
    runTime_.eqns().scalarSources().addSource
    (
        runTime_.value(),
        "t",
        dimTime
    );
    
    //- current time step
    runTime_.eqns().scalarSources().addSource
    (
        runTime_.deltaT().value(),
        "deltaT",
        dimTime
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::admVariableManager::foundStandard(const word& varName) const
{
    forAll(standard_, varIndex)
    {
        if (standard_[varIndex].name() == varName)
        {
            return true;
        }
    }
    return false;
}


bool Foam::admVariableManager::foundImplicit(const word& varName) const
{
    forAll(implicit_, varIndex)
    {
        if (implicit_[varIndex].name() == varName)
        {
            return true;
        }
    }
    return false;
}


bool Foam::admVariableManager::foundDerived(const word& varName) const
{
    forAll(derived_, varIndex)
    {
        if (derived_[varIndex].name() == varName)
        {
            return true;
        }
    }
    return false;
}


bool Foam::admVariableManager::found(const word& varName) const
{
    forAll(all_, varIndex)
    {
        if (all_[varIndex].name() == varName)
        {
            return true;
        }
    }
    return false;
}


const Foam::admStandardVariable& Foam::admVariableManager::lookupStandard
(
    const word& varName
) const
{
    forAll(standard_, varIndex)
    {
        if (standard_[varIndex].name() == varName)
        {
            return standard_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookupStandard")
        << varName << " not found.  Available standard variables are:"
        << tocStandard()
        << abort(FatalError);
    return nullStandard_;
}


Foam::admStandardVariable& Foam::admVariableManager::lookupStandard
(
    const word& varName
)
{
    forAll(standard_, varIndex)
    {
        if (standard_[varIndex].name() == varName)
        {
            return standard_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookupStandard")
        << varName << " not found.  Available standard variables are:"
        << tocStandard()
        << abort(FatalError);
    return nullStandard_;
}


const Foam::admImplicitVariable& Foam::admVariableManager::lookupImplicit
(
    const word& varName
) const
{
    forAll(implicit_, varIndex)
    {
        if (implicit_[varIndex].name() == varName)
        {
            return implicit_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookupImplicit")
        << varName << " not found.  Available implicit variables are:"
        << tocImplicit()
        << abort(FatalError);
    return nullImplicit_;
}


Foam::admImplicitVariable& Foam::admVariableManager::lookupImplicit
(
    const word& varName
)
{
    forAll(implicit_, varIndex)
    {
        if (implicit_[varIndex].name() == varName)
        {
            return implicit_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookupImplicit")
        << varName << " not found.  Available implicit variables are:"
        << tocImplicit()
        << abort(FatalError);
    return nullImplicit_;
}


const Foam::admDerivedVariable& Foam::admVariableManager::lookupDerived
(
    const word& varName
) const
{
    forAll(derived_, varIndex)
    {
        if (derived_[varIndex].name() == varName)
        {
            return derived_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookupDerived")
        << varName << " not found.  Available derived variables are:"
        << tocDerived()
        << abort(FatalError);
    return nullDerived_;
}


Foam::admDerivedVariable& Foam::admVariableManager::lookupDerived
(
    const word& varName
)
{
    forAll(derived_, varIndex)
    {
        if (derived_[varIndex].name() == varName)
        {
            return derived_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookupDerived")
        << varName << " not found.  Available derived variables are:"
        << tocDerived()
        << abort(FatalError);
    return nullDerived_;
}


const Foam::admVariable& Foam::admVariableManager::lookup
(
    const word& varName
) const
{
    forAll(all_, varIndex)
    {
        if (all_[varIndex].name() == varName)
        {
            return all_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookup")
        << varName << " not found.  Available variables are:"
        << toc()
        << abort(FatalError);
    return nullStandard_;
}


Foam::admVariable& Foam::admVariableManager::lookup
(
    const word& varName
)
{
    forAll(all_, varIndex)
    {
        if (all_[varIndex].name() == varName)
        {
            return all_[varIndex];
        }
    }
    FatalErrorIn("admVariableManager::lookup")
        << varName << " not found.  Available variables are:"
        << toc()
        << abort(FatalError);
    return nullStandard_;
}


Foam::wordList Foam::admVariableManager::tocStandard() const
{
    wordList returnMe(standard_.size());
    forAll(standard_, i)
    {
        returnMe[i] = standard_[i].name();
    }
    return returnMe;
}


Foam::wordList Foam::admVariableManager::tocImplicit() const
{
    wordList returnMe(implicit_.size());
    forAll(implicit_, i)
    {
        returnMe[i] = implicit_[i].name();
    }
    return returnMe;
}


Foam::wordList Foam::admVariableManager::tocDerived() const
{
    wordList returnMe(derived_.size());
    forAll(derived_, i)
    {
        returnMe[i] = derived_[i].name();
    }
    return returnMe;
}


Foam::wordList Foam::admVariableManager::toc() const
{
    wordList returnMe(all_.size());
    forAll(all_, i)
    {
        returnMe[i] = all_[i].name();
    }
    return returnMe;
}


void Foam::admVariableManager::initialize()
{
    if (initialized_)
    {
        FatalErrorIn("admVariableManager::initialize")
            << "Object already initialized.  This is considered fatal."
            << abort(FatalError);
    }
    forAll(derived_, varIndex)
    {
        derived_[varIndex].initialize();
    }
    initialized_ = true;
}


void Foam::admVariableManager::massConservingApplyLimits
(
    scalarField& sf,
    const scalarField& weights,
    const scalar& upperLimit,
    const scalar& lowerLimit,
    const word& varName
)
{
    scalarField atLimits(sf.size());
    scalar weightedOffset(0.0);
    scalar weightedSum(0.0);
    while ((max(sf) > upperLimit) || (min(sf) < lowerLimit))
    {
        label limitHits(0);
        forAll(sf, cellIndex)
        {
            if (sf[cellIndex] >= upperLimit)
            {
                limitHits++;
                weightedOffset +=
                    (sf[cellIndex] - upperLimit) * weights[cellIndex];
                atLimits[cellIndex] = upperLimit;
                sf[cellIndex] = 0.0;
            }
            else if (sf[cellIndex] <= lowerLimit)
            {
                limitHits++;
                weightedOffset +=
                    (sf[cellIndex] - lowerLimit) * weights[cellIndex];
                atLimits[cellIndex] = lowerLimit;
                sf[cellIndex] = 0.0;
            }
            else
            {
                weightedSum +=
                    sf[cellIndex] * weights[cellIndex];
                atLimits[cellIndex] = 0.0;
            }
        }
        if (limitHits == sf.size())
        {
            if (sf.size() > 1)
            {
                WarningIn("admVariableManager::massConservingApplyLimits")
                    << varName << " out of limits across full field.  Gained "
                    << "mass = " << weightedOffset << endl;
            }
            sf = atLimits;
            return;
        }
        sf = sf * (1 + weightedOffset / weightedSum) + atLimits;
    }
}


void Foam::admVariableManager::massConservingApplyLimits
(
    admStandardVariable& var
)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    massConservingApplyLimits
    (
        var().internalField(),
        mesh_.V(),
        upperLimit,
        lowerLimit,
        var.name()
    );
}


void Foam::admVariableManager::massConservingApplyLimits
(
    admImplicitVariable& var
)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    massConservingApplyLimits
    (
        var().internalField(),
        mesh_.V(),
        upperLimit,
        lowerLimit,
        var.name()
    );
}


void Foam::admVariableManager::massConservingApplyStandardLimits()
{
    forAll(standard_, varIndex)
    {
        massConservingApplyLimits(standard_[varIndex]);
    }
}


void Foam::admVariableManager::massConservingApplyImplicitLimits()
{
    forAll(implicit_, varIndex)
    {
        massConservingApplyLimits(implicit_[varIndex]);
    }
}


void Foam::admVariableManager::applyLimits(admVariable& var)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    scalarField sf
    (
        max(var.evaluateField(), lowerLimit)
    );
    var.assignField
    (
        min(sf, upperLimit)
    );
}


void Foam::admVariableManager::applyLimits(admStandardVariable& var)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    max(var().internalField(), var().internalField(), lowerLimit);
    min(var().internalField(), var().internalField(), upperLimit);
}


void Foam::admVariableManager::applyLimits(admImplicitVariable& var)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    max(var().internalField(), var().internalField(), lowerLimit);
    min(var().internalField(), var().internalField(), upperLimit);
}


void Foam::admVariableManager::applyLimits
(
    admVariable& var,
    label cellIndex
)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    scalar s
    (
        max(var.evaluate(cellIndex), lowerLimit)
    );
    var.assign
    (
        min(s, upperLimit),
        cellIndex
    );
}


void Foam::admVariableManager::applyLimits
(
    admStandardVariable& var,
    label cellIndex
)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    var().internalField()[cellIndex] = max
    (
        var().internalField()[cellIndex],
        lowerLimit
    );
    var().internalField()[cellIndex] = min
    (
        var().internalField()[cellIndex],
        upperLimit
    );
}


void Foam::admVariableManager::applyLimits
(
    admImplicitVariable& var,
    label cellIndex
)
{
    scalar upperLimit
    (
        var.found("upperLimit")
      ? readScalar(var.lookup("upperLimit"))
      : VGREAT
    );
    scalar lowerLimit
    (
        var.found("lowerLimit")
      ? readScalar(var.lookup("lowerLimit"))
      : scalar(0)
    );
    var().internalField()[cellIndex] = max
    (
        var().internalField()[cellIndex],
        lowerLimit
    );
    var().internalField()[cellIndex] = min
    (
        var().internalField()[cellIndex],
        upperLimit
    );
}


void Foam::admVariableManager::applyStandardLimits()
{
    forAll(standard_, varIndex)
    {
        applyLimits(standard_[varIndex]);
    }
}


void Foam::admVariableManager::applyStandardLimits(label cellIndex)
{
    forAll(standard_, varIndex)
    {
        applyLimits(standard_[varIndex], cellIndex);
    }
}


void Foam::admVariableManager::applyImplicitLimits()
{
    forAll(implicit_, varIndex)
    {
        applyLimits(implicit_[varIndex]);
    }
}


void Foam::admVariableManager::applyImplicitLimits(label cellIndex)
{
    forAll(implicit_, varIndex)
    {
        applyLimits(implicit_[varIndex], cellIndex);
    }
}


void Foam::admVariableManager::applyDerivedLimits()
{
    forAll(derived_, varIndex)
    {
        applyLimits(derived_[varIndex]);
    }
}


void Foam::admVariableManager::applyDerivedLimits(label cellIndex)
{
    forAll(derived_, varIndex)
    {
        applyLimits(derived_[varIndex], cellIndex);
    }
}


void Foam::admVariableManager::saveStandardInternalFields
(
    PtrList<scalarField>& saveSpot
)
{
    if (saveSpot.size() != nStandard())
    {
        saveSpot.clear();
        saveSpot.setSize(nStandard());
        forAll(saveSpot, varIndex)
        {
            saveSpot.set
            (
                varIndex,
                new scalarField
                (
                    standard_[varIndex]().internalField()
                )
            );
        }
    }
    else
    {
        forAll(standard_, varIndex)
        {
            saveSpot[varIndex] = standard_[varIndex]().internalField();
        }
    }
}


void Foam::admVariableManager::saveImplicitInternalFields
(
    PtrList<scalarField>& saveSpot
)
{
    if (saveSpot.size() != nImplicit())
    {
        saveSpot.clear();
        saveSpot.setSize(nImplicit());
        forAll(saveSpot, varIndex)
        {
            saveSpot.set
            (
                varIndex,
                new scalarField
                (
                    implicit_[varIndex]().internalField()
                )
            );
        }
    }
    else
    {
        forAll(implicit_, varIndex)
        {
            saveSpot[varIndex] = implicit_[varIndex]().internalField();
        }
    }
}


void Foam::admVariableManager::saveDerivedInternalFields
(
    PtrList<scalarField>& saveSpot
)
{
    if (saveSpot.size() != nDerived())
    {
        saveSpot.clear();
        saveSpot.setSize(nDerived());
        forAll(saveSpot, varIndex)
        {
            saveSpot.set
            (
                varIndex,
                new scalarField
                (
                    derived_[varIndex].evaluateField()
                )
            );
        }
    }
    else
    {
        forAll(derived_, varIndex)
        {
            saveSpot[varIndex] = derived_[varIndex].evaluateField();
        }
    }
}


void Foam::admVariableManager::saveAllInternalFields
(
    PtrList<scalarField>& saveSpot
)
{
    if (saveSpot.size() != nAll())
    {
        saveSpot.clear();
        saveSpot.setSize(nAll());
        forAll(saveSpot, varIndex)
        {
            saveSpot.set
            (
                varIndex,
                new scalarField
                (
                    all_[varIndex].evaluateField()
                )
            );
        }
    }
    else
    {
        forAll(all_, varIndex)
        {
            saveSpot[varIndex] = all_[varIndex].evaluateField();
        }
    }
}


void Foam::admVariableManager::updateAverages() const
{
    if (model_.outputReactionAverages())
    {
        averagesDict_.clear();
        
        forAll(standard_, varIndex)
        {
            averagesDict_.add
            (
                standard_[varIndex].name(),
                standard_[varIndex]().weightedAverage(mesh_.V()).value()
            );
        }

        forAll(implicit_, varIndex)
        {
            averagesDict_.add
            (
                implicit_[varIndex].name(),
                implicit_[varIndex]().weightedAverage(mesh_.V()).value()
            );
        }

        forAll(derived_, varIndex)
        {
            if (!derived_[varIndex].suppressOutput())
            {
                averagesDict_.add
                (
                    derived_[varIndex].name(),
                    derived_[varIndex].evaluateDimensionedField()
                        .weightedAverage(mesh_.V()).value()
                );
            } // end if not suppressing output
        } // end forAll derived
    } // end if output averages
}


void Foam::admVariableManager::correctStandardBoundaryConditions()
{
    forAll(standard_, varIndex)
    {
        standard_[varIndex].correctBoundaryConditions();
    }
}


void Foam::admVariableManager::correctImplicitBoundaryConditions()
{
    forAll(implicit_, varIndex)
    {
        implicit_[varIndex].correctBoundaryConditions();
    }
}


void Foam::admVariableManager::correctDerivedBoundaryConditions()
{
    forAll(derived_, varIndex)
    {
        derived_[varIndex].correctBoundaryConditions();
    }
}


void Foam::admVariableManager::correctAllBoundaryConditions()
{
    forAll(all_, varIndex)
    {
        all_[varIndex].correctBoundaryConditions();
    }
}


void Foam::admVariableManager::saveState(const label slot)
{
    if ((slot > nStates()) || (slot < 0))
    {
        FatalErrorIn("admVariableManager::saveState")
            << "saveState slot " << slot << " out of range: 0, "
            << nStates() << ". Can only save over existing "
            << "slots or create one at end."
            << abort(FatalError);
    }
    if (slot == nStates())
    {
        standardRestore_.setSize(slot + 1);
        implicitRestore_.setSize(slot + 1);
        standardRestore_.set
        (
            slot,
            new PtrList<volScalarField>()
        );
        implicitRestore_.set
        (
            slot,
            new PtrList<volScalarField>()
        );
    }
    else
    {
        standardRestore_[slot].clear();
        implicitRestore_[slot].clear();
    }
    standardRestore_[slot].setSize(standard_.size());
    implicitRestore_[slot].setSize(implicit_.size());
    forAll(standard_, varIndex)
    {
        standardRestore_[slot].set
        (
            varIndex,
            new volScalarField(standard_[varIndex]())
        );
    }
    forAll(implicit_, varIndex)
    {
        implicitRestore_[slot].set
        (
            varIndex,
            new volScalarField(implicit_[varIndex]())
        );
    }
}


void Foam::admVariableManager::clearState(const label slot)
{
    if ((slot >= nStates()) || (slot < 0))
    {
        WarningIn("admVariableManager::clearState")
            << "Attempting to clear non-existent slot " << slot << ". "
            << "Ignoring. Range = 0, " << nStates() << endl;
        return;
    }
    standardRestore_[slot].clear();
    implicitRestore_[slot].clear();
}


void Foam::admVariableManager::loadState(const label slot)
{
    if ((slot >= nStates()) || (slot < 0))
    {
        FatalErrorIn("admVariableManager::loadState")
            << "Load slot " << slot << " out of range 0, " << nStates()
            << abort(FatalError);
    }
    if (!standardRestore_[slot].size())
    {
        FatalErrorIn("admVariableManager::loadState")
            << "Attempting to load empty save state " << slot
            << abort(FatalError);
    }

    // Loading state is not as straightforward as saving because we have to
    // remap all the .oldTime() references.
    forAll(standard_, varIndex)
    {
        volScalarField * standardPtr(&standard_[varIndex]());
        volScalarField * restorePtr(&standardRestore_[slot][varIndex]);
        label nOldTimes(standard_[varIndex]().nOldTimes());
        for (label delorian(0); delorian < nOldTimes; delorian++)
        {
            *standardPtr = *restorePtr;
            standardPtr = &(standardPtr->oldTime());
            restorePtr = &(restorePtr->oldTime());
        }
        *standardPtr = *restorePtr;
    }
    forAll(implicit_, varIndex)
    {
        volScalarField * implicitPtr(&implicit_[varIndex]());
        volScalarField * restorePtr(&implicitRestore_[slot][varIndex]);
        label nOldTimes(implicit_[varIndex]().nOldTimes());
        for (label delorian(0); delorian < nOldTimes; delorian++)
        {
            *implicitPtr = *restorePtr;
            implicitPtr = &(implicitPtr->oldTime());
            restorePtr = &(restorePtr->oldTime());
        }
        *implicitPtr = *restorePtr;
    }
}


Foam::label Foam::admVariableManager::nStates() const
{
    return standardRestore_.size();
}


bool Foam::admVariableManager::validState(const label slot) const
{
    // Is slot out-of-range?
    if ((slot >= nStates()) || (slot < 0))
    {
        return false;
    }
    
    // Does save slot have no data?
    if (!standardRestore_[slot].size())
    {
        return false;
    }
    
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
