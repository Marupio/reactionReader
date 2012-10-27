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

#include "admTemperatureDependentExponentialRatioCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        admTemperatureDependentExponentialRatioCoefficient,
        0
    );

    addToRunTimeSelectionTable
    (
        admCoefficient,
        admTemperatureDependentExponentialRatioCoefficient,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admTemperatureDependentExponentialRatioCoefficient
    ::admTemperatureDependentExponentialRatioCoefficient
(
    admTime& runTime,
    const fvMesh& mesh,
    admVariableManager& admVars,
    const dictionary& dict,
    const word& coeffName
)
:
    admCoefficient(admVars, coeffName, dimless),
    lastValue_
    (
        scalarField
        (
            mesh_.nCells()
        )
    ),
    lastMilestone_(-1),

    A_
    (
        runTime_.readScalarOrEquation
        (
            dict.subDict(coeffName).lookup("A"),
            word(coeffName + "(A)")
        )
    ),
    B_
    (
        runTime_.readScalarOrEquation
        (
            dict.subDict(coeffName).lookup("B"),
            word(coeffName + "(B)")
        )
    ),
    T_base_
    (
        runTime_.readScalarOrEquation
        (
            dict.subDict(coeffName).lookup("T_base"),
            word(coeffName + "(T_base)")
        )
    ),
    
    varT_
    (
        admVars_.lookup
        (
            word
            (
                dict.subDict(coeffName).lookup("T_var")
            )
        )
    ),
    derivedVarTPtr_
    (
        admVars_.foundDerived
        (
            word
            (
                dict.subDict(coeffName).lookup("T_var")
            )
        )
      ? &admVars_.lookupDerived
        (
            word
            (
                dict.subDict(coeffName).lookup("T_var")
            )
        )
      : NULL
    )
{
    if (dict.subDict(coeffName).found("dimensions"))
    {
        dimensions_.reset(dict.subDict(coeffName).lookup("dimensions"));
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
    Foam::admTemperatureDependentExponentialRatioCoefficient::evaluate
(
    label cellIndex
) const
{
    if (runTime_.milestone() != lastMilestone_)
    {
        return A_ * exp
        (
            B_ * (1 / T_base_ - 1 / varT_.evaluate(cellIndex))
        );
    }
    return lastValue_[cellIndex];
}


Foam::tmp<scalarField>
    Foam::admTemperatureDependentExponentialRatioCoefficient::evaluateField()
    const
{
    if (runTime_.milestone() != lastMilestone_)
    {
        scalarField invTBase(mesh_.nCells(), 1 / T_base_);
        lastValue_ = A_ * exp
        (
            B_ * (invTBase - 1 / varT_.evaluateField())
        );
        lastMilestone_ = runTime_.milestone();
    }

    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            lastValue_
        )
    );
    return tReturnMe;
}


Foam::scalar Foam::admTemperatureDependentExponentialRatioCoefficient::ddt
(
    label cellIndex
) const
{
    if (derivedVarTPtr_)
    {
        const scalar& T(varT_.evaluate(cellIndex));
        return -A_ * B_ / T / T * exp
        (
            B_ * (1 / T_base_ + 1 / T)
        ) * derivedVarTPtr_->ddt(cellIndex);
    }
    return 0.0;
}


Foam::tmp<scalarField>
    Foam::admTemperatureDependentExponentialRatioCoefficient::ddtField() const
{
    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            mesh_.nCells()
        )
    );
    scalarField& returnMe(tReturnMe());

    if (derivedVarTPtr_)
    {
        const scalarField& T(varT_.evaluateField());
        scalarField invTBase(mesh_.nCells(), 1 / T_base_);
        returnMe = -A_ * B_ / T / T * exp
        (
            B_ * (invTBase + 1 / T)
        ) * derivedVarTPtr_->ddtField();
    }
    else
    {
        returnMe = 0.0;
    }
    return tReturnMe;
}


bool Foam::admTemperatureDependentExponentialRatioCoefficient::ddtNonZero()
    const
{
    if (derivedVarTPtr_)
    {
        return (derivedVarTPtr_->ddtNonZero());
    }
    return false;
}


Foam::scalar Foam::admTemperatureDependentExponentialRatioCoefficient::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    if (ddyNonZero(wrtVar))
    {
        const scalar& T(varT_.evaluate(cellIndex));
        return -A_ * B_ / T / T * exp
        (
            B_ * (1 / T_base_ + 1 / T)
        ) * varT_.ddy(wrtVar, cellIndex);
    }
    return 0.0;
}


Foam::tmp<scalarField>
    Foam::admTemperatureDependentExponentialRatioCoefficient::ddyField
(
    const admVariable& wrtVar
) const
{
    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            mesh_.nCells()
        )
    );
    scalarField& returnMe(tReturnMe());
    
    if (ddyNonZero(wrtVar))
    {
        const scalarField& T(varT_.evaluateField());
        scalarField invTBase(mesh_.nCells(), 1 / T_base_);
        returnMe = -A_ * B_ / T / T * exp
        (
            B_ * (invTBase + 1 / T)
        ) * varT_.ddyField(wrtVar);
    }
    else
    {
        returnMe = 0.0;
    }
    return tReturnMe;
}


bool Foam::admTemperatureDependentExponentialRatioCoefficient::ddyNonZero
(
    const admVariable& wrtVar
) const
{
    return varT_.ddyNonZero(wrtVar);
}

// ************************************************************************* //
