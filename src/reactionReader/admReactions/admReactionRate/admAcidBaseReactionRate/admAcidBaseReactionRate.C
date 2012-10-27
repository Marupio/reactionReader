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

#include "admAcidBaseReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admAcidBaseReactionRate, 0);

    addToRunTimeSelectionTable
    (
        admReactionRate,
        admAcidBaseReactionRate,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineAdmJumpTableCoefficientOrVariable(coeff_kAB_, admAcidBaseReactionRate)
defineAdmJumpTableCoefficientOrVariable(coeff_Ka_, admAcidBaseReactionRate)
defineAdmJumpTableCoefficientOrVariable(var_Sm_, admAcidBaseReactionRate)
defineAdmJumpTableCoefficientOrVariable(var_S_, admAcidBaseReactionRate)
defineAdmJumpTableCoefficientOrVariable(var_Shp_, admAcidBaseReactionRate)


admJumpTableDerivative(coeff_kAB_, admAcidBaseReactionRate)
{
    return scalar
    (
        (
            var_Sm_.evaluate(cellIndex) * var_Shp_.evaluate(cellIndex)
          - coeff_Ka_.evaluate(cellIndex) * var_S_.evaluate(cellIndex)
          + (1.0 - using_S_hvar_) * var_Sm_.evaluate(cellIndex)
          * coeff_Ka_.evaluate(cellIndex)
        ) * coeff_kAB_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_kAB_, admAcidBaseReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            (
                var_Sm_.evaluateField() * var_Shp_.evaluateField()
              - coeff_Ka_.evaluateField() * var_S_.evaluateField()
              + (1.0 - using_S_hvar_) * var_Sm_.evaluateField()
              * coeff_Ka_.evaluateField()
            ) * coeff_kAB_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(coeff_Ka_, admAcidBaseReactionRate)
{
    return scalar
    (
        coeff_kAB_.evaluate(cellIndex)
      * (
            (1 - using_S_hvar_) * var_Sm_.evaluate(cellIndex)
          - var_S_.evaluate(cellIndex)
        ) * coeff_Ka_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_Ka_, admAcidBaseReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_kAB_.evaluateField()
          * (
                (1 - using_S_hvar_) * var_Sm_.evaluateField()
              - var_S_.evaluateField()
            ) * coeff_Ka_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_Sm_, admAcidBaseReactionRate)
{
    return scalar
    (
        coeff_kAB_.evaluate(cellIndex)
      * (
            var_Shp_.evaluate(cellIndex)
          + (1 - using_S_hvar_) * coeff_Ka_.evaluate(cellIndex)
        ) * var_Sm_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_Sm_, admAcidBaseReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_kAB_.evaluateField()
          * (
                var_Shp_.evaluateField()
              + (1 - using_S_hvar_) * coeff_Ka_.evaluateField()
            ) * var_Sm_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_S_, admAcidBaseReactionRate)
{
    return scalar
    (
       -coeff_kAB_.evaluate(cellIndex) * coeff_Ka_.evaluate(cellIndex)
      * var_S_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_S_, admAcidBaseReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
           -coeff_kAB_.evaluateField() * coeff_Ka_.evaluateField()
          * var_S_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_Shp_, admAcidBaseReactionRate)
{
    return scalar
    (
        coeff_kAB_.evaluate(cellIndex) * var_Sm_.evaluate(cellIndex)
      * var_Shp_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_Shp_, admAcidBaseReactionRate)
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_kAB_.evaluateField() * var_Sm_.evaluateField()
          * var_Shp_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admAcidBaseReactionRate::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admAcidBaseReactionRate::ddyFieldFns_0
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

Foam::admAcidBaseReactionRate::admAcidBaseReactionRate
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

    initializeAdmReactionRateCoefficient(coeff_kAB_, k_AB),
    initializeAdmReactionRateCoefficient(coeff_Ka_, K_a),

    initializeAdmReactionRateVariable(var_Sm_, S_var_m),
    var_S_
    (
        admVars_.lookup
        (
            word
            (
                using_S_hvar_ == 1.0
              ? dict.subDict(reactionName).subDict("rate").lookup("S_hvar")
              : dict.subDict(reactionName).subDict("rate").lookup("S_var")
            )
        )
    ),
    ddyTablevar_S_(admVars_.nAll()),
    ddyFieldTablevar_S_(admVars_.nAll()),

    initializeAdmReactionRateVariable(var_Shp_, S_H_p),

    using_S_hvar_
    (
        dict.subDict(reactionName).subDict("rate").found("S_hvar")
      ? 1.0 : 0.0
    )
{
    if
    (
        (dimensionSet::debug)
     && (
            (
                var_Sm_.dimensions() * var_Shp_.dimensions()
             != coeff_Ka_.dimensions() * var_S_.dimensions()
            )
         || (coeff_Ka_.dimensions() != var_Shp_.dimensions())
        )
    )
    {
        WarningIn
        (
            "admAcidBaseReactionRate::admAcidBaseReactionRate"
        )
            << "Dimension error thrown in the reaction rate calculation "
            << "for reaction " << reactionName << endl;
    }

    dimensions_.reset
    (
        coeff_kAB_.dimensions()
      * (
            var_Sm_.dimensions() * var_Shp_.dimensions()
          - coeff_Ka_.dimensions() * var_S_.dimensions()
        )
    );

    allocateAdmCoefficientOrVariable(coeff_kAB_, admAcidBaseReactionRate);
    allocateAdmCoefficientOrVariable(coeff_Ka_, admAcidBaseReactionRate);
    allocateAdmCoefficientOrVariable(var_Sm_, admAcidBaseReactionRate);
    allocateAdmCoefficientOrVariable(var_S_, admAcidBaseReactionRate);
    allocateAdmCoefficientOrVariable(var_Shp_, admAcidBaseReactionRate);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admAcidBaseReactionRate::uninhibited
(
    label cellIndex
) const
{
    scalar K_a(coeff_Ka_.evaluate(cellIndex));

    return scalar
    (
        coeff_kAB_.evaluate(cellIndex)
      * (
            var_Sm_.evaluate(cellIndex)
          * (using_S_hvar_ * K_a + var_Shp_.evaluate(cellIndex))
          - K_a * var_S_.evaluate(cellIndex)
        )
    );
}


Foam::tmp<scalarField>
    Foam::admAcidBaseReactionRate::uninhibitedField() const
{
    scalarField K_a(coeff_Ka_.evaluateField());

    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_kAB_.evaluateField()
          * (
                var_Sm_.evaluateField()
              * (using_S_hvar_ * K_a + var_Shp_.evaluateField())
              - K_a * var_S_.evaluateField()
            )
        )
    );
}


Foam::scalar Foam::admAcidBaseReactionRate::uninhibitedDdy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveDerivative(coeff_kAB_)
      + admJumpTableRetrieveDerivative(coeff_Ka_)
      + admJumpTableRetrieveDerivative(var_Sm_)
      + admJumpTableRetrieveDerivative(var_S_)
      + admJumpTableRetrieveDerivative(var_Shp_);      
}


Foam::tmp<scalarField> Foam::admAcidBaseReactionRate::uninhibitedDdyField
(
    const admVariable& wrtVar
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveFieldDerivative(coeff_kAB_)
      + admJumpTableRetrieveFieldDerivative(coeff_Ka_)
      + admJumpTableRetrieveFieldDerivative(var_Sm_)
      + admJumpTableRetrieveFieldDerivative(var_S_)
      + admJumpTableRetrieveFieldDerivative(var_Shp_);      
}


bool Foam::admAcidBaseReactionRate::uninhibitedDdyNonZero
(
    const admVariable& wrtVar
) const
{
    return
    (
        coeff_kAB_.ddyNonZero(wrtVar)
     || coeff_Ka_.ddyNonZero(wrtVar)
     || var_S_.ddyNonZero(wrtVar)
     || var_Sm_.ddyNonZero(wrtVar)
     || var_Shp_.ddyNonZero(wrtVar)
    );
}

// ************************************************************************* //
