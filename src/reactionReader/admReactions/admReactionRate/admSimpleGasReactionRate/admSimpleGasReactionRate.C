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

#include "admSimpleGasReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admSimpleGasReactionRate, 0);

    addToRunTimeSelectionTable
    (
        admReactionRate,
        admSimpleGasReactionRate,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineAdmJumpTableCoefficientOrVariable(coeff_kLa_, admSimpleGasReactionRate)
defineAdmJumpTableCoefficientOrVariable(coeff_n_, admSimpleGasReactionRate)
defineAdmJumpTableCoefficientOrVariable(coeff_KH_, admSimpleGasReactionRate)
defineAdmJumpTableCoefficientOrVariable(var_S_, admSimpleGasReactionRate)
defineAdmJumpTableCoefficientOrVariable(var_P_, admSimpleGasReactionRate)


admJumpTableDerivative(coeff_kLa_, admSimpleGasReactionRate)
{
    return scalar
    (
        (
            var_S_.evaluate(cellIndex)
          - coeff_n_.evaluate(cellIndex) * coeff_KH_.evaluate(cellIndex)
          * var_P_.evaluate(cellIndex)
        ) * coeff_kLa_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_kLa_, admSimpleGasReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            (
                var_S_.evaluateField()
              - coeff_n_.evaluateField() * coeff_KH_.evaluateField()
              * var_P_.evaluateField()
            ) * coeff_kLa_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_n_, admSimpleGasReactionRate)
{
    return scalar
    (
       -coeff_kLa_.evaluate(cellIndex) * coeff_KH_.evaluate(cellIndex)
      * var_P_.evaluate(cellIndex) * coeff_n_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_n_, admSimpleGasReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_kLa_.evaluateField() * coeff_KH_.evaluateField()
          * var_P_.evaluateField() * coeff_n_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_KH_, admSimpleGasReactionRate)
{
    return scalar
    (
       -coeff_kLa_.evaluate(cellIndex) * coeff_n_.evaluate(cellIndex)
      * var_P_.evaluate(cellIndex) * coeff_KH_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_KH_, admSimpleGasReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_kLa_.evaluateField() * coeff_n_.evaluateField()
          * var_P_.evaluateField() * coeff_KH_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_S_, admSimpleGasReactionRate)
{
    return scalar
    (
        coeff_kLa_.evaluate(cellIndex) * var_S_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_S_, admSimpleGasReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_kLa_.evaluateField() * var_S_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_P_, admSimpleGasReactionRate)
{
    return scalar
    (
       -coeff_kLa_.evaluate(cellIndex) * coeff_n_.evaluate(cellIndex)
      * coeff_KH_.evaluate(cellIndex) * var_P_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_P_, admSimpleGasReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_kLa_.evaluateField() * coeff_n_.evaluateField()
          * coeff_KH_.evaluateField() * var_P_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admSimpleGasReactionRate::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admSimpleGasReactionRate::ddyFieldFns_0
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

Foam::admSimpleGasReactionRate::admSimpleGasReactionRate
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

    initializeAdmReactionRateCoefficient(coeff_kLa_, k_L_a),
    initializeAdmReactionRateCoefficient(coeff_n_, K_H),
    initializeAdmReactionRateCoefficient(coeff_KH_, n),

    initializeAdmReactionRateVariable(var_S_, S_var),
    initializeAdmReactionRateVariable(var_P_, p_var)
{
    if
    (
        (dimensionSet::debug)
     && (
            var_S_.dimensions()
         !=
            (
                coeff_n_.dimensions()
              * coeff_KH_.dimensions()
              * var_P_.dimensions()
            )
        )
    )
    {
        WarningIn("admSimpleGasReactionRate::admSimpleGasReactionRate")
            << "Dimension error thrown in the reaction rate calculation "
            << "for reaction " << reactionName << endl;
    }
    dimensions_.reset
    (
        coeff_kLa_.dimensions() *
        (
            var_S_.dimensions()
          + coeff_n_.dimensions()
          * coeff_KH_.dimensions()
          * var_P_.dimensions()
        )
    );
    allocateAdmCoefficientOrVariable(coeff_kLa_, admSimpleGasReactionRate);
    allocateAdmCoefficientOrVariable(coeff_n_, admSimpleGasReactionRate);
    allocateAdmCoefficientOrVariable(coeff_KH_, admSimpleGasReactionRate);
    allocateAdmCoefficientOrVariable(var_S_, admSimpleGasReactionRate);
    allocateAdmCoefficientOrVariable(var_P_, admSimpleGasReactionRate);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admSimpleGasReactionRate::uninhibited
(
    label cellIndex
) const
{
    return scalar
    (
        coeff_kLa_.evaluate(cellIndex)
      * (
            var_S_.evaluate(cellIndex)
          - coeff_n_.evaluate(cellIndex)
          * coeff_KH_.evaluate(cellIndex)
          * var_P_.evaluate(cellIndex)
        )
    );
}


Foam::tmp<scalarField>
    Foam::admSimpleGasReactionRate::uninhibitedField() const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_kLa_.evaluateField()
          * (
                var_S_.evaluateField()
              - coeff_n_.evaluateField()
              * coeff_KH_.evaluateField()
              * var_P_.evaluateField()
            )
        )
    );
}


Foam::scalar Foam::admSimpleGasReactionRate::uninhibitedDdy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveDerivative(coeff_kLa_)
      + admJumpTableRetrieveDerivative(coeff_n_)
      + admJumpTableRetrieveDerivative(coeff_KH_)
      + admJumpTableRetrieveDerivative(var_S_)
      + admJumpTableRetrieveDerivative(var_P_);      
}


Foam::tmp<scalarField> Foam::admSimpleGasReactionRate::uninhibitedDdyField
(
    const admVariable& wrtVar
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveFieldDerivative(coeff_kLa_)
      + admJumpTableRetrieveFieldDerivative(coeff_n_)
      + admJumpTableRetrieveFieldDerivative(coeff_KH_)
      + admJumpTableRetrieveFieldDerivative(var_S_)
      + admJumpTableRetrieveFieldDerivative(var_P_);      
}


bool Foam::admSimpleGasReactionRate::uninhibitedDdyNonZero
(
    const admVariable& wrtVar
) const
{
    return
        coeff_kLa_.ddyNonZero(wrtVar)
     || coeff_n_.ddyNonZero(wrtVar)
     || coeff_KH_.ddyNonZero(wrtVar)
     || var_S_.ddyNonZero(wrtVar)
     || var_P_.ddyNonZero(wrtVar);
}

// ************************************************************************* //
