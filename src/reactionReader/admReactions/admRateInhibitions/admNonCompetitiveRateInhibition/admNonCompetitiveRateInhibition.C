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

#include "admNonCompetitiveRateInhibition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admNonCompetitiveRateInhibition, 0);

    addToRunTimeSelectionTable
    (
        admRateInhibition,
        admNonCompetitiveRateInhibition,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineAdmJumpTableCoefficientOrVariable
(
    coeff_KI_,
    admNonCompetitiveRateInhibition
)
defineAdmJumpTableCoefficientOrVariable
(
    var_SI_,
    admNonCompetitiveRateInhibition
)


admJumpTableDerivative(coeff_KI_, admNonCompetitiveRateInhibition)
{
    scalar denomRoot
    (
        coeff_KI_.evaluate(cellIndex) + var_SI_.evaluate(cellIndex)
    );
    return scalar
    (
        var_SI_.evaluate(cellIndex) / denomRoot / denomRoot
      * coeff_KI_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(coeff_KI_, admNonCompetitiveRateInhibition)
{
    scalarField denomRoot
    (
        coeff_KI_.evaluateField() + var_SI_.evaluateField()
    );
    return tmp<scalarField>
    (
        new scalarField
        (
            var_SI_.evaluateField() / denomRoot / denomRoot
          * coeff_KI_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_SI_, admNonCompetitiveRateInhibition)
{
    scalar KI(coeff_KI_.evaluate(cellIndex));
    scalar denomRoot
    (
        KI + var_SI_.evaluate(cellIndex)
    );
    return scalar
    (
       -KI / denomRoot / denomRoot * var_SI_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_SI_, admNonCompetitiveRateInhibition)
{
    scalarField KI(coeff_KI_.evaluateField());
    scalarField denomRoot
    (
        KI + var_SI_.evaluateField()
    );
    return tmp<scalarField>
    (
        new scalarField
        (
           -KI / denomRoot / denomRoot * var_SI_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admNonCompetitiveRateInhibition::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admNonCompetitiveRateInhibition::ddyFieldFns_0
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

Foam::admNonCompetitiveRateInhibition::
    admNonCompetitiveRateInhibition
(
    admTime& runTime,
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    const dictionary& dict,
    const word& inhibitionName
)
:
    admRateInhibition(admVars, admCoeffs, inhibitionName),
    
    initializeAdmRateInhibitionCoefficient(coeff_KI_, K_I),
    initializeAdmRateInhibitionVariable(var_SI_, S_I)
{
    if
    (
        (dimensionSet::debug)
     && (coeff_KI_.dimensions() != var_SI_.dimensions())
    )
    {
        WarningIn
        (
            "admNonCompetitiveRateInhibition"
            "::admNonCompetitieRateInhibition"
        )
            << "Dimension error thrown for inhibition " << inhibitionName
            << endl;
    }
    dimensionSet testDimensions(coeff_KI_.dimensions());
    testDimensions =  testDimensions + var_SI_.dimensions();
    allocateAdmCoefficientOrVariable
    (
        coeff_KI_,
        admNonCompetitiveRateInhibition
    );
    allocateAdmCoefficientOrVariable
    (
        var_SI_,
        admNonCompetitiveRateInhibition
    );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admNonCompetitiveRateInhibition::inhibition
(
    label cellIndex
) const
{
    scalar K_I(coeff_KI_.evaluate(cellIndex));
    return scalar
    (
        K_I / (K_I + var_SI_.evaluate(cellIndex))
    );
}


Foam::tmp<scalarField>
    Foam::admNonCompetitiveRateInhibition::inhibitionField() const
{
    scalarField K_I(coeff_KI_.evaluateField());
    return tmp<scalarField>
    (
        new scalarField
        (
            K_I / (K_I + var_SI_.evaluateField())
        )
    );
}


bool Foam::admNonCompetitiveRateInhibition::inhibitionNonZero() const
{
    return true;
}


Foam::scalar Foam::admNonCompetitiveRateInhibition::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveDerivative(coeff_KI_)
      + admJumpTableRetrieveDerivative(var_SI_);
}


Foam::tmp<scalarField> Foam::admNonCompetitiveRateInhibition::ddyField
(
    const admVariable& wrtVar
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveFieldDerivative(coeff_KI_)
      + admJumpTableRetrieveFieldDerivative(var_SI_);
}


bool Foam::admNonCompetitiveRateInhibition::ddyNonZero
(
    const admVariable& wrtVar
) const
{
    return
    (
        coeff_KI_.ddyNonZero(wrtVar)
     || var_SI_.ddyNonZero(wrtVar)
    );
}

// ************************************************************************* //
