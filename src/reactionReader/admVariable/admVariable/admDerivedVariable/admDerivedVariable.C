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

#include "admDerivedVariable.H"
#include "admVariableManager.H"
#include "IOstreams.H"
#include "token.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admDerivedVariable, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::admDerivedVariable::createJacobianIndices()
{
    if (!jdict_) return;
    wordList wrtVars(jdict_->toc());
    forAll (wrtVars, i)
    {
        if (!admVars_.found(wrtVars[i]))
        {
            WarningIn("admDerivedVariable::createJacobianIndices")
                << name() << " has a Jacobian function defined "
                << "for unknown variable " << wrtVars[i] << ".  Function "
                << "is being ignored.";
            continue;
        }
        admVariable& var(admVars_.lookup(wrtVars[i]));
        label newIndex(jVars_.size());
        jVars_.setSize(newIndex + 1);
        jNames_.setSize(newIndex + 1);
        jEqns_.setSize(newIndex + 1);

        jVars_.set
        (
            newIndex,
            &var
        );
        word jacobianName
        (
            "d(" + name() + ")|d(" + wrtVars[i] + ")"
        );
        jNames_[newIndex] = jacobianName;

        jEqns_[newIndex] = runTime_.eqns().readEquation
        (
            equation
            (
                jdict_->lookup
                (
                    wrtVars[i]
                ),
                jacobianName
            )
        );
    }
    delete jdict_;
}


void Foam::admDerivedVariable::createReverseFunctionIndices()
{
    if (!rdict_) return;
    wordList wrtVarNames(rdict_->toc());
    wordList order;
    if (rdict_->found("order"))
    {
        order = wordList(rdict_->lookup("order"));
        if (order.size() != (wrtVarNames.size() - 1))
        {
            FatalErrorIn("admDerivedVariable::createReverseFunctionIndices")
                << name() << " reverse function order list size mismatch. "
                << "'order' lists " << order.size() << ", and there are "
                << label(wrtVarNames.size() - 1) << " reverse functions "
                << "defined."
                << abort(FatalError);
        }
    }
    else
    {
        // Assume order doesn't matter
        order = wrtVarNames;
    }
    
    rVars_.setSize(order.size());
    rNames_.setSize(order.size());
    rEqns_.setSize(order.size());
    label offset(0);
    forAll(order, i)
    {
        if (!rdict_->found(order[i]))
        {
            FatalErrorIn("admDerivedVariable::createReverseFunctionIndices")
                << name() << " reverse function order list doesn't match "
                << "defined reverse functions.  order list: " << order
                << ", reverse functions: " << wrtVarNames
                << abort(FatalError);
        }
        if (!admVars_.found(order[i]))
        {
            WarningIn("admDerivedVariable::createReverseFunctionIndices")
                << name() << " has a reverse function defined "
                << "for unknown variable " << order[i] << ".  Function "
                << "is being ignored.";
            offset++;
            continue;
        }
        rVars_.set
        (
            i - offset,
            &admVars_.lookup(order[i])
        );
        word reverseName
        (
            "reverse(" + name() + ")solvedFor(" + order[i] + ")"
        );
        rNames_[i - offset] = reverseName;
        rEqns_[i - offset] = runTime_.eqns().readEquation
        (
            equation
            (
                rdict_->lookup
                (
                    order[i]
                ),
                reverseName
            )
        );
    }
    if (offset)
    {
        label newSize(rVars_.size() - offset);
        rVars_.setSize(newSize);
        rNames_.setSize(newSize);
        rEqns_.setSize(newSize);
    }
    delete rdict_;
}


void Foam::admDerivedVariable::findChainedVariables
(
    const label& searchDepth,
    UPtrList<admVariable>& varList,
    labelListList& eqnIndices,
    labelList& operations
) const
{
    // Create index lists for this variable - will be empty unless it is a
    // derived variable
    varList = jVars_;
    eqnIndices.setSize(jEqns_.size());
    forAll(eqnIndices, i)
    {
        eqnIndices[i].setSize(1);
        eqnIndices[i][0] = jEqns_[i];
    }

    operations.setSize(varList.size());
    forAll(operations, i)
    {
        operations[i] =
            runTime_.eqns()[eqnIndices[i][0]].size();
    }

    if (searchDepth > 1)
    {
        forAll(jVars_, i)
        {
            UPtrList<admVariable> nextVars;
            labelListList nextEqnIndices;
            labelList nextOperations;
            const admVariable& nextVariable(jVars_[i]);
            if (nextVariable.varType() == admVariable::vtderived)
            {
                const admDerivedVariable& nextDerivedVariable
                (
                    admVars_.derived(nextVariable.localIndex())
                );
                nextDerivedVariable.findChainedVariables
                (
                    searchDepth - 1,
                    nextVars,
                    nextEqnIndices,
                    nextOperations
                );
            }
            // add operations
            forAll(nextOperations, j)
            {
                nextOperations[j] += operations[i];
                label newEqnIndex(nextEqnIndices[j].size());
                nextEqnIndices[j].setSize(newEqnIndex + 1);
                nextEqnIndices[j][newEqnIndex] = eqnIndices[i][0];
            }
            appendChainedLists
            (
                nextVars,
                nextEqnIndices,
                nextOperations,
                varList,
                eqnIndices,
                operations
            );
        }
    }
}


void Foam::admDerivedVariable::appendChainedLists
(
    UPtrList<admVariable> nextVars,
    labelListList& nextEqnIndices,
    labelList& nextOperations,
    UPtrList<admVariable> currentVars,
    labelListList& currentEqnIndices,
    labelList& currentOperations
) const
{
    // Go through all that needs to be appended
    forAll(nextVars, appendMe)
    {
        bool duplicateFound(false);
        
        // Check for duplicates
        forAll(currentVars, testMe)
        {
            if (nextVars[appendMe] == currentVars[testMe])
            {
                // Duplicate detected, overwrite if fewer operations
                duplicateFound = true;
                if (nextOperations[appendMe] < currentOperations[testMe])
                {
                    currentEqnIndices[testMe] = nextEqnIndices[appendMe];
                    currentOperations[testMe] = nextOperations[appendMe];
                }
                break;
            }
        }
        if (!duplicateFound)
        {
            label newIndex(currentVars.size());
            currentVars.setSize(newIndex + 1);
            currentEqnIndices.setSize(newIndex + 1);
            currentOperations.setSize(newIndex + 1);
            
            currentVars.set(newIndex, &nextVars[appendMe]);
            currentEqnIndices[newIndex] = nextEqnIndices[appendMe];
            currentOperations[newIndex] = nextOperations[appendMe];
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admDerivedVariable::admDerivedVariable
(
    admTime& runTime,
    admVariableManager& admVars
)
:
    admVariable
    (
        runTime,
        admVars,
        dictionary(),
        -1,
        -1
    ),
    dataPtr_(NULL)
{}


Foam::admDerivedVariable::admDerivedVariable
(
    admTime& runTime,
    admVariableManager& admVars,
    const dictionary& dict,
    const label globalIndex,
    const label localIndex,
    const word& name
)
:
    admVariable
    (
        runTime,
        admVars,
        dict,
        globalIndex,
        localIndex
    ),

    suppressOutput_
    (
        dict.found("suppressOutput")
      ? dict.lookup("suppressOutput")
      : Switch(false)
    ),

    dataPtr_
    (
        new volScalarField
        (
            IOobject
            (
                name,
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                suppressOutput_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(name, dimless, 0.0)
        )
    ),

    dimensionsKnown_(false),
    
    // read main function equation, store equation index
    funEqn_
    (
        runTime_.eqns().readEquation
        (
            equation(lookup("function"), name)
        )
    ),

    // read ddt equation if present, store equation index, or -1 if not present
    ddtEqn_
    (
        found("ddt")
      ? runTime_.eqns().readEquation
        (
            equation(lookup("ddt"), ddtName())
        )
      : label(-1)
    ),
    
    // The reverseFunctions can't be read in until admVars knows all the
    // variables.  For now we store a pointer to the dictionary.
    rdict_
    (
        found("reverseFunctions")
      ? new dictionary(subDict("reverseFunctions"))
      : NULL
    ),
    rVars_(0),
    rNames_(0),
    rEqns_(0),

    // The Jacobian can't be read in until admVars knows all the variables
    // For now we store a pointer to the dictionary
    jdict_
    (
        found("Jacobian")
      ? new dictionary(subDict("Jacobian"))
      : NULL
    ),
    jVars_(0),
    jNames_(0),
    jEqns_(0),
    cVars_(0),
    cEqns_(0),
    lastMilestone_(-1)
{
    equationVariable& eqVar(* this);
    runTime_.eqns().addSource
    (
        eqVar
    );

    varType_ = admVariable::vtderived;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::admDerivedVariable::name() const
{
    return dataPtr_().name();
}

const Foam::dimensionSet& Foam::admDerivedVariable::dimensions() const
{
    if (!dimensionsKnown_)
    {
        FatalErrorIn("admDerivedVariable::dimensions")
            << "Dimensions still unknown for variable " << dataPtr_().name()
            << "; must initialize the variable first."
            << abort(FatalError);
    }
    return dataPtr_().dimensions();
}


void Foam::admDerivedVariable::correctBoundaryConditions()
{
    dataPtr_().correctBoundaryConditions();
}


void Foam::admDerivedVariable::createIndices()
{
    createJacobianIndices();
    createReverseFunctionIndices();
}


void Foam::admDerivedVariable::initialize()
{
    dataPtr_().dimensions().reset
    (
        runTime_.eqns().evaluateDimensions(funEqn_)
    );
    dimensionsKnown_ = true;
}


void Foam::admDerivedVariable::createChainRuleIndices
(
    const label chainRuleSearchDepth
)
{
    // Dummy list
    labelList operations;

    findChainedVariables
    (
        chainRuleSearchDepth,
        cVars_,
        cEqns_,
        operations
    );
}


void Foam::admDerivedVariable::assign
(
    const scalar& scalarIn,
    label cellIndex
)
{
    dataPtr_().internalField()[cellIndex] = scalarIn;
    forAll(rEqns_, varIndex)
    {
        rVars_[varIndex].assign
        (
            runTime_.eqns().evaluateScalar(rEqns_[varIndex], cellIndex),
            cellIndex
        );
    }
}


void Foam::admDerivedVariable::assignDims
(
    const dimensionSet& dimsIn
)
{
    dimensions() = dimsIn;
    forAll(rEqns_, varIndex)
    {
        rVars_[varIndex].assignDims
        (
            runTime_.eqns().evaluateDimensions(rEqns_[varIndex])
        );
    }
}


void Foam::admDerivedVariable::assignDimensioned
(
    const dimensionedScalar& dScalarIn,
    label cellIndex
)
{
    assign(dScalarIn.value(), cellIndex);
    assignDims(dScalarIn.dimensions());
}


void Foam::admDerivedVariable::assignField
(
    const scalarField& sfIn
)
{
    dataPtr_().internalField() = sfIn;
    forAll(rEqns_, varIndex)
    {
        scalarField result(sfIn.size());
        runTime_.eqns().evaluateScalarField(result, rEqns_[varIndex]);
        rVars_[varIndex].assignField(result);
    }
}


void Foam::admDerivedVariable::assignDimensionedField
(
    const dimensionedScalarField& dsfIn
)
{
    assignDims(dsfIn.dimensions());
    assignField(dsfIn.field());
}


const Foam::scalar& Foam::admDerivedVariable::evaluate
(
    label cellIndex
) const
{
    if (lastMilestone_ != runTime_.milestone())
    {
        dataPtr_().internalField()[cellIndex] =
            runTime_.eqns().evaluateScalar(funEqn_, cellIndex);
    }
    return dataPtr_().internalField()[cellIndex];
    // do not set new milestone - only for full field calculations.
}


Foam::word Foam::admDerivedVariable::evaluateName() const
{
    return name();
}


const Foam::dimensionSet& Foam::admDerivedVariable::evaluateDims() const
{
    return dimensions();
}


Foam::dimensionedScalar Foam::admDerivedVariable::evaluateDimensioned
(
    label cellIndex
) const
{
    return dimensionedScalar
    (
        evaluateName(),
        evaluateDims(),
        evaluate(cellIndex)
    );
}


const Foam::scalarField& Foam::admDerivedVariable::evaluateField() const
{
    scalarField& varField(dataPtr_().internalField());
    if (lastMilestone_ != runTime_.milestone())
    {
        runTime_.eqns().evaluateScalarField(varField, funEqn_);
        lastMilestone_ = runTime_.milestone();
    }
    return varField;
}


const Foam::dimensionedScalarField&
    Foam::admDerivedVariable::evaluateDimensionedField() const
{
    evaluateDims();
    evaluateField();
    return dataPtr_().dimensionedInternalField();
}


const Foam::volScalarField&
    Foam::admDerivedVariable::evaluateGeometricField() const
{
    evaluateDims();
    evaluateField();
    return dataPtr_();
}


bool Foam::admDerivedVariable::evaluateNonZero() const
{
    return true;
}


Foam::scalar Foam::admDerivedVariable::ddt(label cellIndex) const
{
    if (ddtEqn_ >= 0)
    {
        return runTime_.eqns().evaluateScalar(ddtEqn_, cellIndex);
    }
    return 0;
}


Foam::word Foam::admDerivedVariable::ddtName() const
{
    return word("ddt(" + name() + ")");
}


Foam::dimensionSet Foam::admDerivedVariable::ddtDims() const
{
    return dimensionSet(dimensions() / dimTime);
}


Foam::dimensionedScalar Foam::admDerivedVariable::ddtDimensioned
(
    label cellIndex
) const
{
    return dimensionedScalar
    (
        ddtName(),
        ddtDims(),
        ddt(cellIndex)
    );
}


Foam::tmp<scalarField> Foam::admDerivedVariable::ddtField() const
{
    tmp<scalarField> tReturnMe
    (
        new scalarField(admVars_.mesh().nCells())
    );
    scalarField& returnMe(tReturnMe());
    if (ddtEqn_ >= 0)
    {
        runTime_.eqns().evaluateScalarField(returnMe, ddtEqn_);
    }
    else
    {
        returnMe = 0.0;
    }
    return tReturnMe;
}


Foam::tmp<dimensionedScalarField>
    Foam::admDerivedVariable::ddtDimensionedField() const
{
    return tmp<dimensionedScalarField>
    (
        new dimensionedScalarField
        (
            IOobject
            (
                ddtName(),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            ddtDims(),
            ddtField()
        )
    );
}


bool Foam::admDerivedVariable::ddtNonZero() const
{
    return (ddtEqn_ >= 0);
}


Foam::scalar Foam::admDerivedVariable::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Check itself
    if (wrtVar == *this)
    {
        return scalar(1.0);
    }

    // Check its Jacobian lists
    forAll(jVars_, varIndex)
    {
        if (jVars_[varIndex] == wrtVar)
        {
            return runTime_.eqns().evaluateScalar
            (
                jEqns_[varIndex],
                cellIndex
            );
        }
    }
    
    // Check its chainRule lists
    forAll(cVars_, i)
    {
        if (cVars_[i] == wrtVar)
        {
            scalar returnMe(1);
            forAll(cEqns_[i], j)
            {
                returnMe *= runTime_.eqns().evaluateScalar
                (
                    cEqns_[i][j],
                    cellIndex
                );
            }
            return returnMe;
        }
    }
    return scalar(0.0);
}


Foam::word Foam::admDerivedVariable::ddyName
(
    const admVariable& wrtVar
) const
{
    return word("d(" + name() + ")|d(" + wrtVar.name() + ")");
}


Foam::dimensionSet Foam::admDerivedVariable::ddyDims
(
    const admVariable& wrtVar
) const
{
    return dimensions() / wrtVar.dimensions();
}


Foam::dimensionedScalar Foam::admDerivedVariable::ddyDimensioned
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return dimensionedScalar
    (
        ddyName(wrtVar),
        ddyDims(wrtVar),
        ddy(wrtVar, cellIndex)
    );
}


Foam::tmp<scalarField> Foam::admDerivedVariable::ddyField
(
    const admVariable& wrtVar
) const
{

    // Check itself
    if (wrtVar == *this)
    {
        return tmp<scalarField>
        (
            new scalarField
            (
                mesh_.nCells(),
                1.0
            )
        );
    }

    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            mesh_.nCells()
        )
    );
    scalarField& returnMe(tReturnMe());

    // Check Jacobian lists
    forAll(jVars_, varIndex)
    {
        if (jVars_[varIndex] == wrtVar)
        {
            runTime_.eqns().evaluateScalarField(returnMe, jEqns_[varIndex]);
            return tReturnMe;
        }
    }

    // Check chain rule lists
    forAll(cVars_, i)
    {
        if (cVars_[i] == wrtVar)
        {
            returnMe = 1.0;

            forAll(cEqns_[i], j)
            {
                scalarField result(returnMe.size());
                runTime_.eqns().evaluateScalarField(result, cEqns_[i][j]);
                returnMe *= result;
            }
            return tReturnMe;
        }
    }
    returnMe = 0.0;
    return tReturnMe;
}


Foam::tmp<dimensionedScalarField> Foam::admDerivedVariable::ddyDimensionedField
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
                ddyName(wrtVar),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            ddyDims(wrtVar),
            ddyField(wrtVar)
        )
    );
}


bool Foam::admDerivedVariable::ddyNonZero(const admVariable& wrtVar) const
{
    // Check itself
    if (wrtVar == *this)
    {
        return true;
    }
    
    // Check its Jacobian lists
    forAll(jVars_, i)
    {
        if (jVars_[i] == wrtVar)
        {
            return true;
        }
    }

    // Check its chainRule lists
    forAll(cVars_, i)
    {
        if (cVars_[i] == wrtVar)
        {
            return true;
        }
    }
    return false;
}


const volScalarField& admDerivedVariable::operator()() const
{
    return dataPtr_();
}


Foam::label Foam::admDerivedVariable::lookupComponentIndex(const word) const
{
    return label(0);
}


Foam::scalar Foam::admDerivedVariable::evaluateScalar
(
    const label componentIndex,
    label cellIndex,
    const label geoIndex
) const
{
    return evaluate(cellIndex);
}


void Foam::admDerivedVariable::evaluateScalarField
(
    scalarField& result,
    const label componentIndex,
    const label geoIndex
) const
{
    result = evaluateField();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
