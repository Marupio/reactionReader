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

#include "admStandardVariable.H"
#include "admVariableManager.H"
#include "IOstreams.H"
#include "token.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admStandardVariable, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admStandardVariable::admStandardVariable
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


Foam::admStandardVariable::admStandardVariable
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

    changedByUdf_
    (
        dict.found("changedByUdf")
      ? dict.lookup("changedByUdf")
      : Switch(false)
    )
{
    runTime_.eqns().scalarSources().addSource
    (
        dataPtr_()
    );

    varType_ = admVariable::vtstandard;
    
    if (convergence_ < 0)
    {
        convergence_ = admVars_.defaultStandardConvergence();
        if (convergence_ < 0)
        {
            FatalIOErrorIn("admStandardVariable::admStandardVariable", dict)
                << "A positive convergence criterion must be defined here or "
                << "in defaults."
                << exit(FatalIOError);
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::admStandardVariable::name() const
{
    return dataPtr_().name();
}


const Foam::dimensionSet& Foam::admStandardVariable::dimensions() const
{
    return dataPtr_().dimensions();
}


void Foam::admStandardVariable::assign
(
    const scalar& scalarIn,
    label cellIndex
)
{
    dataPtr_().internalField()[cellIndex] = scalarIn;
}


void Foam::admStandardVariable::assignDims
(
    const dimensionSet& dimsIn
)
{
    dataPtr_().dimensions() = dimsIn;
}


void Foam::admStandardVariable::assignDimensioned
(
    const dimensionedScalar& dScalarIn,
    label cellIndex
)
{
    assignDims(dScalarIn.dimensions());
    assign(dScalarIn.value(), cellIndex);
}


void Foam::admStandardVariable::assignField
(
    const scalarField& sfIn
)
{
    dataPtr_().internalField() = sfIn;
}


void Foam::admStandardVariable::assignDimensionedField
(
    const dimensionedScalarField& dsfIn
)
{
    assignDims(dsfIn.dimensions());
    assignField(dsfIn.field());
}


const Foam::scalar& Foam::admStandardVariable::evaluate(label cellIndex) const
{
    return dataPtr_().internalField()[cellIndex];
}


Foam::word Foam::admStandardVariable::evaluateName() const
{
    return name();
}


const Foam::dimensionSet& Foam::admStandardVariable::evaluateDims() const
{
    return dimensions();
}


Foam::dimensionedScalar Foam::admStandardVariable::evaluateDimensioned
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


const scalarField& Foam::admStandardVariable::evaluateField() const
{
    return dataPtr_().internalField();
}


const Foam::dimensionedScalarField&
    Foam::admStandardVariable::evaluateDimensionedField() const
{
    return dataPtr_().dimensionedInternalField();
}


const Foam::volScalarField&
    Foam::admStandardVariable::evaluateGeometricField() const
{
    return dataPtr_();
}


bool Foam::admStandardVariable::evaluateNonZero() const
{
    return true;
}


Foam::scalar Foam::admStandardVariable::ddy
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

Foam::word Foam::admStandardVariable::ddyName(const admVariable& wrtVar) const
{
    return word("d" + name() + "|d" + wrtVar.name());
}


Foam::dimensionSet Foam::admStandardVariable::ddyDims
(
    const admVariable& wrtVar
) const
{
    return dimensionSet
    (
        dimensions() / wrtVar.dimensions()
    );
}

Foam::dimensionedScalar Foam::admStandardVariable::ddyDimensioned
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

Foam::tmp<scalarField> Foam::admStandardVariable::ddyField
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

Foam::tmp<dimensionedScalarField> Foam::admStandardVariable::ddyDimensionedField
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

bool Foam::admStandardVariable::ddyNonZero(const admVariable& wrtVar) const
{
    // Check itself
    if (wrtVar == *this)
    {
        return true;
    }
    return false;
}


void Foam::admStandardVariable::correctBoundaryConditions()
{
    dataPtr_().correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
