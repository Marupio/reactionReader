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

#include "admImplicitVariable.H"
#include "admVariableManager.H"
#include "IOstreams.H"
#include "token.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admImplicitVariable, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admImplicitVariable::admImplicitVariable
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
    dataPtr_(NULL),
    gamma_
    (
        "null(gamma)",
        dimless,
        0.0
    )
{}


Foam::admImplicitVariable::admImplicitVariable
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

    dataPtr_
    (
        new volScalarField
        (
            IOobject
            (
                name,
                runTime_.timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),

    gamma_
    (
        lookup("diffusion")
    ),

    convergence_
    (
        found("convergence")
      ? readScalar(lookup("convergence"))
      : -VGREAT
    ),

    minErrorScale_
    (
        found("minErrorScale")
      ? readScalar(lookup("minErrorScale"))
      : SMALL
    ),

    autoSolve_
    (
        found("autoSolve")
      ? lookup("autoSolve")
      : Switch(false)
    ),
    autoSolveConvergence_
    (
        found("autoSolveConvergence")
      ? readScalar(lookup("autoSolveConvergence"))
      : convergence_
    ),
    autoSolveMaxIter_
    (
        found("autoSolveMaxIter")
      ? readLabel(lookup("autoSolveMaxIter"))
      : label(10000)
    )
{
    runTime_.eqns().scalarSources().addSource
    (
        dataPtr_()
    );

    varType_ = admVariable::vtimplicit;

    if (convergence_ < 0)
    {
        convergence_ = admVars_.defaultImplicitConvergence();
        if (convergence_ < 0)
        {
            FatalIOErrorIn("admImplicitVariable::admImplicitVariable", dict)
                << "A positive convergence criterion must be defined here or "
                << "in defaults."
                << exit(FatalIOError);
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::admImplicitVariable::name() const
{
    return dataPtr_().name();
}


const Foam::dimensionSet& Foam::admImplicitVariable::dimensions() const
{
    return dataPtr_().dimensions();
}


void Foam::admImplicitVariable::assign
(
    const scalar& scalarIn,
    label cellIndex
)
{
    dataPtr_().internalField()[cellIndex] = scalarIn;
}


void Foam::admImplicitVariable::assignDims
(
    const dimensionSet& dimsIn
)
{
    dataPtr_().dimensions() = dimsIn;
}


void Foam::admImplicitVariable::assignDimensioned
(
    const dimensionedScalar& dScalarIn,
    label cellIndex
)
{
    assignDims(dScalarIn.dimensions());
    assign(dScalarIn.value(), cellIndex);
}


void Foam::admImplicitVariable::assignField
(
    const scalarField& sfIn
)
{
    dataPtr_().internalField() = sfIn;
}


void Foam::admImplicitVariable::assignDimensionedField
(
    const dimensionedScalarField& dsfIn
)
{
    assignDims(dsfIn.dimensions());
    assignField(dsfIn.field());
}


const Foam::scalar& Foam::admImplicitVariable::evaluate(label cellIndex) const
{
    return dataPtr_().internalField()[cellIndex];
}


Foam::word Foam::admImplicitVariable::evaluateName() const
{
    return name();
}


const Foam::dimensionSet& Foam::admImplicitVariable::evaluateDims() const
{
    return dimensions();
}


Foam::dimensionedScalar Foam::admImplicitVariable::evaluateDimensioned
(
    label cellIndex
) const
{
    return dimensionedScalar
    (
        name(),
        evaluateDims(),
        evaluate(cellIndex)
    );
}


const scalarField& Foam::admImplicitVariable::evaluateField() const
{
    return dataPtr_().internalField();
}


const Foam::dimensionedScalarField&
    Foam::admImplicitVariable::evaluateDimensionedField() const
{
    return dataPtr_().dimensionedInternalField();
}


const Foam::volScalarField&
    Foam::admImplicitVariable::evaluateGeometricField() const
{
    return dataPtr_();
}


bool Foam::admImplicitVariable::evaluateNonZero() const
{
    return true;
}


Foam::scalar Foam::admImplicitVariable::ddy
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
    else
    {
        return scalar(0.0);
    }
}

Foam::word Foam::admImplicitVariable::ddyName(const admVariable& wrtVar) const
{
    return word("d" + name() + "|d" + wrtVar.name());
}


Foam::dimensionSet Foam::admImplicitVariable::ddyDims
(
    const admVariable& wrtVar
) const
{
    return dimensionSet
    (
        dimensions() / wrtVar.dimensions()
    );
}

Foam::dimensionedScalar Foam::admImplicitVariable::ddyDimensioned
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

Foam::tmp<scalarField> Foam::admImplicitVariable::ddyField
(
    const admVariable& wrtVar
) const
{
    tmp<scalarField> tReturnMe
    (
        new scalarField(mesh_.nCells())
    );
    scalarField& returnMe(tReturnMe());

    if (wrtVar == *this)
    {
        returnMe = 1.0;
    }
    else
    {
        returnMe = 0.0;
    }

    return tReturnMe;
}

Foam::tmp<dimensionedScalarField> Foam::admImplicitVariable::ddyDimensionedField
(
    const admVariable& wrtVar
) const
{
    tmp<dimensionedScalarField> tReturnMe
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
    return tReturnMe;
}

bool Foam::admImplicitVariable::ddyNonZero(const admVariable& wrtVar) const
{
    // Check itself
    if (wrtVar == *this)
    {
        return true;
    }
    return false;
}


void Foam::admImplicitVariable::correctBoundaryConditions()
{
    dataPtr_().correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
