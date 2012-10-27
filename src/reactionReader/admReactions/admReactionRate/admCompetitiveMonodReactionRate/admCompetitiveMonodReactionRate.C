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

#include "admCompetitiveMonodReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admCompetitiveMonodReactionRate, 0);

    addToRunTimeSelectionTable
    (
        admReactionRate,
        admCompetitiveMonodReactionRate,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineAdmJumpTableCoefficientOrVariable
(
    coeff_km_,
    admCompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    coeff_Ks_,
    admCompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    coeff_KI_,
    admCompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    var_S_,
    admCompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    var_X_,
    admCompetitiveMonodReactionRate
)
defineAdmJumpTableCoefficientOrVariable
(
    var_SI_,
    admCompetitiveMonodReactionRate
)


admJumpTableDerivative(coeff_km_, admCompetitiveMonodReactionRate)
{
    scalar KI(coeff_KI_.evaluate(cellIndex));
    return scalar
    (
        KI * var_S_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      / (
            coeff_Ks_.evaluate(cellIndex)
          * (KI + var_SI_.evaluate(cellIndex))
          + var_S_.evaluate(cellIndex) * KI
        ) * coeff_km_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_km_, admCompetitiveMonodReactionRate)
{
    scalarField KI(coeff_KI_.evaluateField());
    return tmp<scalarField>
    (
        new scalarField
        (
            KI * var_S_.evaluateField() * var_X_.evaluateField()
          / (
                coeff_Ks_.evaluateField()
              * (KI + var_SI_.evaluateField())
              + var_S_.evaluateField() * KI
            ) * coeff_km_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_Ks_, admCompetitiveMonodReactionRate)
{
    scalar KI(coeff_KI_.evaluate(cellIndex));
    scalar denomRoot
    (
        coeff_Ks_.evaluate(cellIndex) * (KI + var_SI_.evaluate(cellIndex))
      + var_S_.evaluate(cellIndex) * KI
    );
    return scalar
    (
       -coeff_km_.evaluate(cellIndex) * KI * var_S_.evaluate(cellIndex)
      * var_X_.evaluate(cellIndex) * (KI + var_SI_.evaluate(cellIndex))
      / denomRoot / denomRoot * coeff_Ks_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_Ks_, admCompetitiveMonodReactionRate)
{
    scalarField KI(coeff_KI_.evaluateField());
    scalarField denomRoot
    (
        coeff_Ks_.evaluateField() * (KI + var_SI_.evaluateField())
      + var_S_.evaluateField() * KI
    );
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_km_.evaluateField() * KI * var_S_.evaluateField()
          * var_X_.evaluateField() * (KI + var_SI_.evaluateField())
          / denomRoot / denomRoot * coeff_Ks_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_KI_, admCompetitiveMonodReactionRate)
{
    scalar Ks(coeff_Ks_.evaluate(cellIndex));
    scalar KI(coeff_KI_.evaluate(cellIndex));
    scalar denomRoot
    (
        Ks * (KI + var_SI_.evaluate(cellIndex))
      + var_S_.evaluate(cellIndex) * KI
    );
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * var_S_.evaluate(cellIndex)
      * var_X_.evaluate(cellIndex)
      * (
            (
                Ks * (KI + var_SI_.evaluate(cellIndex))
              + var_S_.evaluate(cellIndex) * KI
            )
          - (
                Ks + var_S_.evaluate(cellIndex)
            ) * KI
        )
      / denomRoot / denomRoot * coeff_KI_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_KI_, admCompetitiveMonodReactionRate)
{
    scalarField Ks(coeff_Ks_.evaluateField());
    scalarField KI(coeff_KI_.evaluateField());
    scalarField denomRoot
    (
        Ks * (KI + var_SI_.evaluateField())
      + var_S_.evaluateField() * KI
    );
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * var_S_.evaluateField()
          * var_X_.evaluateField()
          * (
                (
                    Ks * (KI + var_SI_.evaluateField())
                  + var_S_.evaluateField() * KI
                )
              - (
                    Ks + var_S_.evaluateField()
                ) * KI
            )
          / denomRoot / denomRoot * coeff_KI_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_S_, admCompetitiveMonodReactionRate)
{
    scalar Ks(coeff_Ks_.evaluate(cellIndex));
    scalar KI(coeff_KI_.evaluate(cellIndex));
    scalar denomRoot
    (
        Ks * (KI + var_SI_.evaluate(cellIndex))
      + var_S_.evaluate(cellIndex) * KI
    );
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * KI * var_X_.evaluate(cellIndex)
      * Ks * (KI + var_SI_.evaluate(cellIndex)) / denomRoot / denomRoot
      * var_S_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_S_, admCompetitiveMonodReactionRate)
{
    scalarField Ks(coeff_Ks_.evaluateField());
    scalarField KI(coeff_KI_.evaluateField());
    scalarField denomRoot
    (
        Ks * (KI + var_SI_.evaluateField())
      + var_S_.evaluateField() * KI
    );
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * KI * var_X_.evaluateField()
          * Ks * (KI + var_SI_.evaluateField()) / denomRoot / denomRoot
          * var_S_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_X_, admCompetitiveMonodReactionRate)
{
    scalar KI(coeff_KI_.evaluate(cellIndex));
    return scalar
    (
        coeff_km_.evaluate(cellIndex) * KI
      * var_S_.evaluate(cellIndex) * var_X_.evaluate(cellIndex)
      / (
            coeff_Ks_.evaluate(cellIndex)
          * (KI + var_SI_.evaluate(cellIndex))
          + var_S_.evaluate(cellIndex) * KI
        ) * var_X_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_X_, admCompetitiveMonodReactionRate)
{
    scalarField KI(coeff_KI_.evaluateField());
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * KI
          * var_S_.evaluateField() * var_X_.evaluateField()
          / (
                coeff_Ks_.evaluateField()
              * (KI + var_SI_.evaluateField())
              + var_S_.evaluateField() * KI
            ) * var_X_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_SI_, admCompetitiveMonodReactionRate)
{
    scalar Ks(coeff_Ks_.evaluate(cellIndex));
    scalar KI(coeff_KI_.evaluate(cellIndex));
    scalar denomRoot
    (
        Ks * (KI + var_SI_.evaluate(cellIndex))
      + var_S_.evaluate(cellIndex) * KI
    );
    return scalar
    (
       -coeff_km_.evaluate(cellIndex) * KI * var_S_.evaluate(cellIndex)
      * var_X_.evaluate(cellIndex) * Ks / denomRoot / denomRoot
      * var_SI_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_SI_, admCompetitiveMonodReactionRate)
{
    scalarField Ks(coeff_Ks_.evaluateField());
    scalarField KI(coeff_KI_.evaluateField());
    scalarField denomRoot
    (
        Ks * (KI + var_SI_.evaluateField())
      + var_S_.evaluateField() * KI
    );
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_km_.evaluateField() * KI * var_S_.evaluateField()
          * var_X_.evaluateField() * Ks / denomRoot / denomRoot
          * var_SI_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admCompetitiveMonodReactionRate::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admCompetitiveMonodReactionRate::ddyFieldFns_0
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

Foam::admCompetitiveMonodReactionRate::admCompetitiveMonodReactionRate
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
            (coeff_KI_.dimensions() != var_SI_.dimensions())
         || (coeff_Ks_.dimensions() != var_S_.dimensions())
        )
    )
    {
        WarningIn
        (
            "admCompetitiveMonodReactionRate::admCompetitiveMonodReactionRate"
        )
            << "Dimension error thrown in the reaction rate calculation "
            << "for reaction " << reactionName << endl;
    }

    dimensions_.reset
    (
        coeff_km_.dimensions() * coeff_KI_.dimensions()
      * var_S_.dimensions() * var_X_.dimensions()
      / (
            coeff_Ks_.dimensions()
          * (coeff_KI_.dimensions() + var_SI_.dimensions())
          + var_S_.dimensions() * coeff_KI_.dimensions()
        )
    );
    allocateAdmCoefficientOrVariable
    (
        coeff_km_,
        admCompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        coeff_Ks_,
        admCompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        coeff_KI_,
        admCompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        var_S_,
        admCompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        var_X_,
        admCompetitiveMonodReactionRate
    );
    allocateAdmCoefficientOrVariable
    (
        var_SI_,
        admCompetitiveMonodReactionRate
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admCompetitiveMonodReactionRate::uninhibited
(
    label cellIndex
) const
{
    scalar K_I(coeff_KI_.evaluate(cellIndex));

    return scalar
    (
        coeff_km_.evaluate(cellIndex) * K_I
      * var_S_.evaluate(cellIndex) * var_X_.evaluate(cellIndex) 
      / (
            coeff_Ks_.evaluate(cellIndex)
          * (K_I + var_SI_.evaluate(cellIndex))
          + var_S_.evaluate(cellIndex) * K_I
        )
    );
}


Foam::tmp<scalarField>
    Foam::admCompetitiveMonodReactionRate::uninhibitedField() const
{
    scalarField K_I(coeff_KI_.evaluateField());

    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_km_.evaluateField() * K_I
          * var_S_.evaluateField() * var_X_.evaluateField() 
          / (
                coeff_Ks_.evaluateField()
              * (K_I + var_SI_.evaluateField())
              + var_S_.evaluateField() * K_I
            )
        )
    );
}


Foam::scalar Foam::admCompetitiveMonodReactionRate::uninhibitedDdy
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


Foam::tmp<scalarField>
    Foam::admCompetitiveMonodReactionRate::uninhibitedDdyField
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


bool Foam::admCompetitiveMonodReactionRate::uninhibitedDdyNonZero
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
