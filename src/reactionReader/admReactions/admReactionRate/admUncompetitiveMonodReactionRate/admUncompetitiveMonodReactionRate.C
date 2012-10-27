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

#include "admUncompetitiveMonodReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admUncompetitiveMonodReactionRate, 0);

    addToRunTimeSelectionTable
    (
        admReactionRate,
        admUncompetitiveMonodReactionRate,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineAdmJumpTableCoefficientOrVariable
(
    coeff_km_,
    admUncompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    coeff_Ks_,
    admUncompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    coeff_KI_,
    admUncompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    var_S_,
    admUncompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    var_X_,
    admUncompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    var_SI_,
    admUncompetitiveMonodReactionRate
)


admJumpTableDerivative(coeff_km_, admUncompetitiveMonodReactionRate)
{
    return scalar
    (
        var_X_.evaluate(cellIndex) * var_S_.evaluate(cellIndex)
      * var_SI_.evaluate(cellIndex)
      / (
            coeff_Ks_.evaluate(cellIndex) * var_SI_.evaluate(cellIndex)
          + var_S_.evaluate(cellIndex)
          * (var_SI_.evaluate(cellIndex) + coeff_KI_.evaluate(cellIndex))
        ) * coeff_km_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_km_, admUncompetitiveMonodReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            var_X_.evaluateField() * var_S_.evaluateField()
          * var_SI_.evaluateField()
          / (
                coeff_Ks_.evaluateField() * var_SI_.evaluateField()
              + var_S_.evaluateField()
              * (var_SI_.evaluateField() + coeff_KI_.evaluateField())
            ) * coeff_km_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_Ks_, admUncompetitiveMonodReactionRate)
{
    scalar denomRoot
    (
        coeff_Ks_.evaluate(cellIndex) * var_SI_.evaluate(cellIndex)
      + var_S_.evaluate(cellIndex)
      * (var_SI_.evaluate(cellIndex) + coeff_KI_.evaluate(cellIndex))
    );
    return scalar
    (
       -coeff_km_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      * var_S_.evaluate(cellIndex) * pow(var_SI_.evaluate(cellIndex), 2.0)
      / denomRoot / denomRoot * coeff_Ks_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_Ks_, admUncompetitiveMonodReactionRate)
{
    scalarField denomRoot
    (
        coeff_Ks_.evaluateField() * var_SI_.evaluateField()
      + var_S_.evaluateField()
      * (var_SI_.evaluateField() + coeff_KI_.evaluateField())
    );
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_km_.evaluateField() * var_X_.evaluateField()
          * var_S_.evaluateField() * pow(var_SI_.evaluateField(), 2.0)
          / denomRoot / denomRoot * coeff_Ks_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_KI_, admUncompetitiveMonodReactionRate)
{
    scalar denomRoot
    (
        coeff_Ks_.evaluate(cellIndex) * var_SI_.evaluate(cellIndex)
      + var_S_.evaluate(cellIndex)
      * (var_SI_.evaluate(cellIndex) + coeff_KI_.evaluate(cellIndex))
    );
    return scalar
    (
       -coeff_km_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      * pow(var_S_.evaluate(cellIndex), 2.0) * var_SI_.evaluate(cellIndex)
      / denomRoot / denomRoot * coeff_Ks_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_KI_, admUncompetitiveMonodReactionRate)
{
    scalarField denomRoot
    (
        coeff_Ks_.evaluateField() * var_SI_.evaluateField()
      + var_S_.evaluateField()
      * (var_SI_.evaluateField() + coeff_KI_.evaluateField())
    );
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_km_.evaluateField() * var_X_.evaluateField()
          * pow(var_S_.evaluateField(), 2.0) * var_SI_.evaluateField()
          / denomRoot / denomRoot * coeff_Ks_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_S_, admUncompetitiveMonodReactionRate)
{
    scalar Ks(coeff_Ks_.evaluate(cellIndex));
    scalar denomRoot
    (
        Ks * var_SI_.evaluate(cellIndex)
      + var_S_.evaluate(cellIndex)
      * (var_SI_.evaluate(cellIndex) + coeff_KI_.evaluate(cellIndex))
    );
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      * Ks * pow(var_SI_.evaluate(cellIndex), 2.0)
      / denomRoot / denomRoot * var_S_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_S_, admUncompetitiveMonodReactionRate)
{
    scalarField Ks(coeff_Ks_.evaluateField());
    scalarField denomRoot
    (
        Ks * var_SI_.evaluateField()
      + var_S_.evaluateField()
      * (var_SI_.evaluateField() + coeff_KI_.evaluateField())
    );
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * var_X_.evaluateField()
          * Ks * pow(var_SI_.evaluateField(), 2.0)
          / denomRoot / denomRoot * var_S_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_X_, admUncompetitiveMonodReactionRate)
{
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * var_S_.evaluate(cellIndex)
      * var_SI_.evaluate(cellIndex)
      / (
            coeff_Ks_.evaluate(cellIndex) * var_SI_.evaluate(cellIndex)
          + var_S_.evaluate(cellIndex)
          * (var_SI_.evaluate(cellIndex) + coeff_KI_.evaluate(cellIndex))
        ) * var_X_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_X_, admUncompetitiveMonodReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * var_S_.evaluateField()
          * var_SI_.evaluateField()
          / (
                coeff_Ks_.evaluateField() * var_SI_.evaluateField()
              + var_S_.evaluateField()
              * (var_SI_.evaluateField() + coeff_KI_.evaluateField())
            ) * var_X_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_SI_, admUncompetitiveMonodReactionRate)
{
    scalar KI(coeff_KI_.evaluate(cellIndex));
    scalar denomRoot
    (
        coeff_Ks_.evaluate(cellIndex) * var_SI_.evaluate(cellIndex)
      + var_S_.evaluate(cellIndex) * (var_SI_.evaluate(cellIndex) + KI)
    );
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      * KI * pow(var_S_.evaluate(cellIndex), 2.0)
      / denomRoot / denomRoot * var_SI_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_SI_, admUncompetitiveMonodReactionRate)
{
    scalarField KI(coeff_KI_.evaluateField());
    scalarField denomRoot
    (
        coeff_Ks_.evaluateField() * var_SI_.evaluateField()
      + var_S_.evaluateField() * (var_SI_.evaluateField() + KI)
    );
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * var_X_.evaluateField()
          * KI * pow(var_S_.evaluateField(), 2.0)
          / denomRoot / denomRoot * var_SI_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admUncompetitiveMonodReactionRate::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admUncompetitiveMonodReactionRate::ddyFieldFns_0
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

Foam::admUncompetitiveMonodReactionRate::admUncompetitiveMonodReactionRate
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
    initializeAdmReactionRateCoefficient(coeff_KI_, K_I),
    initializeAdmReactionRateVariable(var_S_, S_var),
    initializeAdmReactionRateVariable(var_X_, X_var),
    initializeAdmReactionRateVariable(var_SI_, S_I)
{
    if
    (
        (dimensionSet::debug)
     && (
            (
                var_SI_.dimensions()
             != coeff_KI_.dimensions()
            )
         || (
                (coeff_Ks_.dimensions() * var_SI_.dimensions())
             != (var_S_.dimensions() * var_SI_.dimensions())
            )
        )
    )
    {
        WarningIn
        (
            "admUncompetitiveMonodReactionRate::"
            "admUncompetitiveMonodReactionRate"
        )
            << "Dimension error thrown in the reaction rate calculation "
            << "for reaction " << reactionName << endl;
    }

    dimensions_.reset
    (
        coeff_km_.dimensions() * var_X_.dimensions()
      * var_S_.dimensions() * var_SI_.dimensions()
      / (
            coeff_Ks_.dimensions() * var_SI_.dimensions() + var_S_.dimensions()
          * (var_SI_.dimensions() + coeff_KI_.dimensions())
        )
    );
    allocateAdmCoefficientOrVariable
    (
        coeff_km_,
        admUncompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        coeff_Ks_,
        admUncompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        coeff_KI_,
        admUncompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        var_S_,
        admUncompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        var_X_,
        admUncompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        var_SI_,
        admUncompetitiveMonodReactionRate
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admUncompetitiveMonodReactionRate::uninhibited
(
    label cellIndex
) const
{
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      * var_S_.evaluate(cellIndex) * var_SI_.evaluate(cellIndex)
      / (
            coeff_Ks_.evaluate(cellIndex) * var_SI_.evaluate(cellIndex)
          + var_S_.evaluate(cellIndex) *
            (var_SI_.evaluate(cellIndex) + coeff_KI_.evaluate(cellIndex))
        )
    );
}


Foam::tmp<scalarField>
    Foam::admUncompetitiveMonodReactionRate::uninhibitedField() const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * var_X_.evaluateField()
          * var_S_.evaluateField() * var_SI_.evaluateField()
          / (
                coeff_Ks_.evaluateField() * var_SI_.evaluateField()
              + var_S_.evaluateField() *
                (var_SI_.evaluateField() + coeff_KI_.evaluateField())
            )
        )
    );
}


Foam::scalar Foam::admUncompetitiveMonodReactionRate::uninhibitedDdy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveDerivative(coeff_km_)
      + admJumpTableRetrieveDerivative(coeff_Ks_)
      + admJumpTableRetrieveDerivative(coeff_KI_)
      + admJumpTableRetrieveDerivative(var_S_)
      + admJumpTableRetrieveDerivative(var_X_)
      + admJumpTableRetrieveDerivative(var_SI_);
}


Foam::tmp<scalarField> Foam::admUncompetitiveMonodReactionRate::uninhibitedDdyField
(
    const admVariable& wrtVar
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveFieldDerivative(coeff_km_)
      + admJumpTableRetrieveFieldDerivative(coeff_Ks_)
      + admJumpTableRetrieveFieldDerivative(coeff_KI_)
      + admJumpTableRetrieveFieldDerivative(var_S_)
      + admJumpTableRetrieveFieldDerivative(var_X_)
      + admJumpTableRetrieveFieldDerivative(var_SI_);
}


bool Foam::admUncompetitiveMonodReactionRate::uninhibitedDdyNonZero
(
    const admVariable& wrtVar
) const
{
    return
    (
        coeff_km_.ddyNonZero(wrtVar)
     || coeff_Ks_.ddyNonZero(wrtVar)
     || coeff_KI_.ddyNonZero(wrtVar)
     || var_S_.ddyNonZero(wrtVar)
     || var_SI_.ddyNonZero(wrtVar)
     || var_X_.ddyNonZero(wrtVar)
    );
}

// ************************************************************************* //
