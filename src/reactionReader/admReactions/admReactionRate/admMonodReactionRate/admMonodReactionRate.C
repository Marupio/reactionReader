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

#include "admMonodReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admMonodReactionRate, 0);

    addToRunTimeSelectionTable
    (
        admReactionRate,
        admMonodReactionRate,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineAdmJumpTableCoefficientOrVariable(coeff_km_, admMonodReactionRate)
defineAdmJumpTableCoefficientOrVariable(coeff_Ks_, admMonodReactionRate)
defineAdmJumpTableCoefficientOrVariable(var_S_, admMonodReactionRate)
defineAdmJumpTableCoefficientOrVariable(var_X_, admMonodReactionRate)


admJumpTableDerivative(coeff_km_, admMonodReactionRate)
{
    return scalar
    (
        var_S_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      / (coeff_Ks_.evaluate(cellIndex) + var_S_.evaluate(cellIndex))
      * coeff_km_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_km_, admMonodReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            var_S_.evaluateField() * var_X_.evaluateField()
          / (coeff_Ks_.evaluateField() + var_S_.evaluateField())
          * coeff_km_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_Ks_, admMonodReactionRate)
{
    scalar denomRoot
    (
        coeff_Ks_.evaluate(cellIndex) + var_S_.evaluate(cellIndex)
    );
    return scalar
    (
       -coeff_km_.evaluate(cellIndex) * var_S_.evaluate(cellIndex)
      * var_X_.evaluate(cellIndex) / denomRoot / denomRoot
      * coeff_Ks_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_Ks_, admMonodReactionRate)
{
    scalarField denomRoot
    (
        coeff_Ks_.evaluateField() + var_S_.evaluateField()
    );
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_km_.evaluateField() * var_S_.evaluateField()
          * var_X_.evaluateField() / denomRoot / denomRoot
          * coeff_Ks_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_S_, admMonodReactionRate)
{
    scalar Ks(coeff_Ks_.evaluate(cellIndex));
    scalar denomRoot
    (
        Ks + var_S_.evaluate(cellIndex)
    );
    return scalar
    (
        Ks * coeff_km_.evaluate(cellIndex)
      * var_X_.evaluate(cellIndex) / denomRoot / denomRoot
      * var_S_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_S_, admMonodReactionRate)
{
    scalarField Ks(coeff_Ks_.evaluateField());
    scalarField denomRoot
    (
        Ks + var_S_.evaluateField()
    );
    return tmp<scalarField>
    (
        new scalarField
        (
            Ks * coeff_km_.evaluateField()
          * var_X_.evaluateField() / denomRoot / denomRoot
          * var_S_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_X_, admMonodReactionRate)
{
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * var_S_.evaluate(cellIndex)
      / (coeff_Ks_.evaluate(cellIndex) + var_S_.evaluate(cellIndex))
      * var_X_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_X_, admMonodReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * var_S_.evaluateField()
          / (coeff_Ks_.evaluateField() + var_S_.evaluateField())
          * var_X_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admMonodReactionRate::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admMonodReactionRate::ddyFieldFns_0
(
    const admVariable& wrtVar
) const
{
    return tmp<scalarField>
    (
        new scalarField(mesh_.nCells(), 0.0)
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admMonodReactionRate::admMonodReactionRate
(
    admTime& runTime,
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    admReactionManager& admReacs,
    const dictionary& dict,
    const word& reactionName
)
:
    admReactionRate(admVars, admCoeffs, admReacs, dict, reactionName),

    initializeAdmReactionRateCoefficient(coeff_km_, k_m),
    initializeAdmReactionRateCoefficient(coeff_Ks_, K_s),
    initializeAdmReactionRateVariable(var_S_, S_var),
    initializeAdmReactionRateVariable(var_X_, X_var)
{
    if
    (
        (dimensionSet::debug)
     && (var_S_.dimensions() != coeff_Ks_.dimensions())
    )
    {
        WarningIn("admMonodReactionRate::admMonodReactionRate")
            << "Dimension error thrown in the reaction rate calculation "
            << "for reaction " << reactionName << endl;
    }

    dimensions_.reset
    (
        coeff_km_.dimensions() * var_S_.dimensions()
        * var_X_.dimensions()
        / (coeff_Ks_.dimensions() + var_S_.dimensions())
    );

    allocateAdmCoefficientOrVariable(coeff_km_, admMonodReactionRate);
    allocateAdmCoefficientOrVariable(coeff_Ks_, admMonodReactionRate);
    allocateAdmCoefficientOrVariable(var_S_, admMonodReactionRate);
    allocateAdmCoefficientOrVariable(var_X_, admMonodReactionRate);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admMonodReactionRate::uninhibited
(
    label cellIndex
) const
{
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * var_S_.evaluate(cellIndex)
      * var_X_.evaluate(cellIndex)
      / (coeff_Ks_.evaluate(cellIndex) + var_S_.evaluate(cellIndex))
    );
}


Foam::tmp<scalarField>
    Foam::admMonodReactionRate::uninhibitedField() const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * var_S_.evaluateField()
          * var_X_.evaluateField()
          / (coeff_Ks_.evaluateField() + var_S_.evaluateField())
        )
    );
}


Foam::scalar Foam::admMonodReactionRate::uninhibitedDdy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveDerivative(coeff_km_)
      + admJumpTableRetrieveDerivative(coeff_Ks_)
      + admJumpTableRetrieveDerivative(var_S_)
      + admJumpTableRetrieveDerivative(var_X_);
}


Foam::tmp<scalarField> Foam::admMonodReactionRate::uninhibitedDdyField
(
    const admVariable& wrtVar
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveFieldDerivative(coeff_km_)
      + admJumpTableRetrieveFieldDerivative(coeff_Ks_)
      + admJumpTableRetrieveFieldDerivative(var_S_)
      + admJumpTableRetrieveFieldDerivative(var_X_);
}


bool Foam::admMonodReactionRate::uninhibitedDdyNonZero
(
    const admVariable& wrtVar
) const
{
    return 
    (
        coeff_km_.ddyNonZero(wrtVar)
     || coeff_Ks_.ddyNonZero(wrtVar)
     || var_S_.ddyNonZero(wrtVar)
     || var_X_.ddyNonZero(wrtVar)
    );
}

// ************************************************************************* //
