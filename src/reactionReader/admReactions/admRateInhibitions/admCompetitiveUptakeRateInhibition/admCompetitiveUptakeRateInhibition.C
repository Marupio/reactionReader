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

#include "admCompetitiveUptakeRateInhibition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admCompetitiveUptakeRateInhibition, 0);

    addToRunTimeSelectionTable
    (
        admRateInhibition,
        admCompetitiveUptakeRateInhibition,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

defineAdmJumpTableCoefficientOrVariable
(
    var_S_,
    admCompetitiveUptakeRateInhibition
)
defineAdmJumpTableCoefficientOrVariable
(
    var_SI_,
    admCompetitiveUptakeRateInhibition
)

admJumpTableDerivative(var_S_, admCompetitiveUptakeRateInhibition)
{
    scalar denomRoot
    (
        stabilise
        (
            var_SI_.evaluate(cellIndex) + var_S_.evaluate(cellIndex),
            stabiliseEpsilon_
        )
    );
    return scalar
    (
        var_SI_.evaluate(cellIndex) / denomRoot / denomRoot
      * var_S_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_S_, admCompetitiveUptakeRateInhibition)
{
    scalarField denomRoot
    (
        stabilise
        (
            var_SI_.evaluateField() + var_S_.evaluateField(),
            stabiliseEpsilon_
        )
    );
    return tmp<scalarField>
    (
        new scalarField
        (
            var_SI_.evaluateField() / denomRoot / denomRoot
          * var_S_.ddyField(wrtVar)
        )
    );
}


admJumpTableDerivative(var_SI_, admCompetitiveUptakeRateInhibition)
{
    scalar denomRoot
    (
        stabilise
        (
            var_SI_.evaluate(cellIndex) + var_S_.evaluate(cellIndex),
            stabiliseEpsilon_
        )
    );
    return scalar
    (
       -var_S_.evaluate(cellIndex) / denomRoot / denomRoot
      * var_SI_.ddy(wrtVar, cellIndex)
    );
}


admJumpTableFieldDerivative(var_SI_, admCompetitiveUptakeRateInhibition)
{
    scalarField denomRoot
    (
        stabilise
        (
            var_SI_.evaluateField() + var_S_.evaluateField(),
            stabiliseEpsilon_
        )
    );
    return tmp<scalarField>
    (
        new scalarField
        (
           -var_S_.evaluateField() / denomRoot / denomRoot
          * var_SI_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admCompetitiveUptakeRateInhibition::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admCompetitiveUptakeRateInhibition::ddyFieldFns_0
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

Foam::admCompetitiveUptakeRateInhibition::
    admCompetitiveUptakeRateInhibition
(
    admTime& runTime,
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    const dictionary& dict,
    const word& inhibitionName
)
:
    admRateInhibition(admVars, admCoeffs, inhibitionName),

    initializeAdmRateInhibitionVariable(var_S_, S_var),
    initializeAdmRateInhibitionVariable(var_SI_, S_I),
    stabiliseEpsilon_
    (
        dict.subDict(inhibitionName).found("stabiliseEpsilon")
      ? readScalar(dict.subDict(inhibitionName).lookup("stabiliseEpsilon"))
      : SMALL
    )
{
    if
    (
        (dimensionSet::debug)
     && (var_S_.dimensions() != var_SI_.dimensions())
    )
    {
        WarningIn
        (
            "admCompetitiveRateInhibition"
            "::admCompetitieRateInhibition"
        )
            << "Dimension error thrown for inhibition " << inhibitionName
            << endl;
    }
    dimensionSet testDimensions(var_S_.dimensions());
    testDimensions =  testDimensions + var_SI_.dimensions();
    allocateAdmCoefficientOrVariable
    (
        var_S_,
        admCompetitiveUptakeRateInhibition
    );
    allocateAdmCoefficientOrVariable
    (
        var_SI_,
        admCompetitiveUptakeRateInhibition
    );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admCompetitiveUptakeRateInhibition::inhibition
(
    label cellIndex
) const
{
    return scalar
    (
        var_S_.evaluate(cellIndex)
      / stabilise
        (
            var_S_.evaluate(cellIndex) + var_SI_.evaluate(cellIndex),
            stabiliseEpsilon_
        )
    );
}


Foam::tmp<scalarField>
    Foam::admCompetitiveUptakeRateInhibition::inhibitionField() const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            var_S_.evaluateField()
          / stabilise
            (
                var_S_.evaluateField() + var_SI_.evaluateField(),
                stabiliseEpsilon_
            )
        )
    );
}


bool Foam::admCompetitiveUptakeRateInhibition::inhibitionNonZero() const
{
    return true;
}


Foam::scalar Foam::admCompetitiveUptakeRateInhibition::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveDerivative(var_S_)
      + admJumpTableRetrieveDerivative(var_SI_);
}


Foam::tmp<scalarField> Foam::admCompetitiveUptakeRateInhibition::ddyField
(
    const admVariable& wrtVar
) const
{
    // Call functions from jump table
    return
        admJumpTableRetrieveFieldDerivative(var_S_)
      + admJumpTableRetrieveFieldDerivative(var_SI_);
}


bool Foam::admCompetitiveUptakeRateInhibition::ddyNonZero
(
    const admVariable& wrtVar
) const
{
    return
    (
        var_S_.ddyNonZero(wrtVar)
     || var_SI_.ddyNonZero(wrtVar)
    );
}

// ************************************************************************* //
