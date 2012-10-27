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

#include "admFirstOrderReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admFirstOrderReactionRate, 0);

    addToRunTimeSelectionTable
    (
        admReactionRate,
        admFirstOrderReactionRate,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::admFirstOrderReactionRate::ddyFns_init_k
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    label globalIndex(wrtVar.globalIndex());
    
    if (coeff_k_.ddyNonZero(wrtVar))
    {
        ddyTable_k_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFns_k;
        ddyFieldTable_k_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFieldFns_k;
    }
    else
    {
        ddyTable_k_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFns_0;
        ddyFieldTable_k_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFieldFns_0;
    }

    // Call the originally requested function
    return (*this.*ddyTable_k_[globalIndex])(wrtVar, cellIndex);
}


Foam::scalar Foam::admFirstOrderReactionRate::ddyFns_init_var
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    label globalIndex(wrtVar.globalIndex());
    
    if (var_.ddyNonZero(wrtVar))
    {
        ddyTable_var_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFns_var;
        ddyFieldTable_var_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFieldFns_var;
    }
    else
    {
        ddyTable_var_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFns_0;
        ddyFieldTable_var_[globalIndex] =
            &Foam::admFirstOrderReactionRate::ddyFieldFns_0;
    }

    // Call the originally requested function
    return (*this.*ddyTable_var_[globalIndex])(wrtVar, cellIndex);
}


Foam::tmp<scalarField> Foam::admFirstOrderReactionRate::ddyFieldFns_init_k
(
    const admVariable& wrtVar
) const
{
    ddyFns_init_k(wrtVar, 0);

    // Call the originally requested function
    return (*this.*ddyFieldTable_k_[wrtVar.globalIndex()])(wrtVar);
}


Foam::tmp<scalarField> Foam::admFirstOrderReactionRate::ddyFieldFns_init_var
(
    const admVariable& wrtVar
) const
{
    ddyFns_init_var(wrtVar, 0);

    // Call the originally requested function
    return (*this.*ddyFieldTable_var_[wrtVar.globalIndex()])(wrtVar);
}


Foam::scalar Foam::admFirstOrderReactionRate::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField> Foam::admFirstOrderReactionRate::ddyFieldFns_0
(
    const admVariable& wrtVar
) const
{
    return tmp<scalarField>
    (
        new scalarField(mesh_.nCells(), 0.0)
    );
}


Foam::scalar Foam::admFirstOrderReactionRate::ddyFns_k
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(var_.evaluate(cellIndex) * coeff_k_.ddy(wrtVar, cellIndex));
}


Foam::tmp<scalarField> Foam::admFirstOrderReactionRate::ddyFieldFns_k
(
    const admVariable& wrtVar
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            var_.evaluateField() * coeff_k_.ddyField(wrtVar)
        )
    );
}


Foam::scalar Foam::admFirstOrderReactionRate::ddyFns_var
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(coeff_k_.evaluate(cellIndex) * var_.ddy(wrtVar, cellIndex));
}


Foam::tmp<scalarField> Foam::admFirstOrderReactionRate::ddyFieldFns_var
(
    const admVariable& wrtVar
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_k_.evaluateField() * var_.ddyField(wrtVar)
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admFirstOrderReactionRate::admFirstOrderReactionRate
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

    coeff_k_
    (
        admCoeffs_.readCoefficient
        (
            dict.subDict(reactionName).subDict("rate"),
            word("k"),
            word(reactionName + "(rate(k))")
        )
    ),
    
    var_
    (
        admVars_.lookup
        (
            word
            (
                dict.subDict(reactionName).subDict("rate").lookup("var")
            )
        )
    ),

    ddyTable_k_(admVars_.nAll()),
    ddyTable_var_(admVars_.nAll()),
    ddyFieldTable_k_(admVars_.nAll()),
    ddyFieldTable_var_(admVars_.nAll())
{
    dimensions_.reset
    (
        coeff_k_.dimensions() * var_.dimensions()
    );

    forAll(ddyTable_k_, index)
    {
        ddyTable_k_.set
        (
            index,
            new ddyFn(&Foam::admFirstOrderReactionRate::ddyFns_init_k)
        );
        ddyTable_var_.set
        (
            index,
            new ddyFn(&Foam::admFirstOrderReactionRate::ddyFns_init_var)
        );
        ddyFieldTable_k_.set
        (
            index,
            new ddyFieldFn
            (
                &Foam::admFirstOrderReactionRate::ddyFieldFns_init_k
            )
        );
        ddyFieldTable_var_.set
        (
            index,
            new ddyFieldFn
            (
                &Foam::admFirstOrderReactionRate::ddyFieldFns_init_var
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admFirstOrderReactionRate::uninhibited
(
    label cellIndex
) const
{
    return scalar
    (
        coeff_k_.evaluate(cellIndex) * var_.evaluate(cellIndex)
    );
}


Foam::tmp<scalarField>
    Foam::admFirstOrderReactionRate::uninhibitedField() const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            coeff_k_.evaluateField() * var_.evaluateField()
        )
    );
}


Foam::scalar Foam::admFirstOrderReactionRate::uninhibitedDdy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call functions from jump table
    return
        (*this.*ddyTable_k_[wrtVar.globalIndex()])(wrtVar, cellIndex)
      + (*this.*ddyTable_var_[wrtVar.globalIndex()])(wrtVar, cellIndex);
}


Foam::tmp<scalarField> Foam::admFirstOrderReactionRate::uninhibitedDdyField
(
    const admVariable& wrtVar
) const
{
    // Call functions from jump table
    return
        (*this.*ddyFieldTable_k_[wrtVar.globalIndex()])(wrtVar)
      + (*this.*ddyFieldTable_var_[wrtVar.globalIndex()])(wrtVar);
}


bool Foam::admFirstOrderReactionRate::uninhibitedDdyNonZero
(
    const admVariable& wrtVar
) const
{
    return var_.ddyNonZero(wrtVar) || coeff_k_.ddyNonZero(wrtVar);
}

// ************************************************************************* //
