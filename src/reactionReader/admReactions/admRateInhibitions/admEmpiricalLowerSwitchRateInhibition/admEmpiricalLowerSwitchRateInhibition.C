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

#include "admEmpiricalLowerSwitchRateInhibition.H"
#include "addToRunTimeSelectionTable.H"
#include "versionSpecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(admEmpiricalLowerSwitchRateInhibition, 0);

    addToRunTimeSelectionTable
    (
        admRateInhibition,
        admEmpiricalLowerSwitchRateInhibition,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::admEmpiricalLowerSwitchRateInhibition::ddyFns_Init
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    label globalIndex(wrtVar.globalIndex());
    if (varPH_.ddyNonZero(wrtVar))
    {
        ddyFns_[globalIndex] =
            &Foam::admEmpiricalLowerSwitchRateInhibition::ddyFns_VarPH;
        ddyFieldFns_[globalIndex] =
            &Foam::admEmpiricalLowerSwitchRateInhibition::ddyFieldFns_VarPH;
    }
    else
    {
        ddyFns_[globalIndex] =
            &Foam::admEmpiricalLowerSwitchRateInhibition::ddyFns_0;
        ddyFieldFns_[globalIndex] =
            &Foam::admEmpiricalLowerSwitchRateInhibition::ddyFieldFns_0;
    }

    // Call the originally requested function
    return (*this.*ddyFns_[globalIndex])(wrtVar, cellIndex);
}


Foam::tmp<scalarField>
    Foam::admEmpiricalLowerSwitchRateInhibition::ddyFieldFns_Init
(
    const admVariable& wrtVar
) const
{
    ddyFns_Init(wrtVar, 0);

    // Call the originally requested function
    return (*this.*ddyFieldFns_[wrtVar.globalIndex()])(wrtVar);
}


Foam::scalar Foam::admEmpiricalLowerSwitchRateInhibition::ddyFns_0
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return scalar(0.0);
}


Foam::tmp<scalarField>
    Foam::admEmpiricalLowerSwitchRateInhibition::ddyFieldFns_0
(
    const admVariable& wrtVar
) const
{
    return tmp<scalarField>
    (
        new scalarField(mesh_.nCells(), 0.0)
    );
}


Foam::scalar Foam::admEmpiricalLowerSwitchRateInhibition::ddyFns_VarPH
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    scalar returnMe(0.0);
    scalar pH_LL(coeff_pH_LL_.evaluate(cellIndex));
    scalar pH_UL(coeff_pH_UL_.evaluate(cellIndex));
    scalar pH_var(varPH_.evaluate(cellIndex));
    scalar multiplier;
    if (S_Hp_)
    {
        multiplier = -1000.0 * log10(MathConstantScope::e) / pH_var;
        pH_var = - log10(pH_var / 1000.0);
    }
    else
    {
        multiplier = 1.0;
    }
    if (pH_var < pH_UL)
    {
        scalar denominatorRootRecip
        (
            1.0 / (pH_UL - pH_LL)
        );
        returnMe =
            6.0 * pH_UL * (pH_var - pH_UL)
          * denominatorRootRecip * denominatorRootRecip
          * exp
            (
                -3.0 * sqr
                (
                    (pH_var - pH_UL) * denominatorRootRecip
                )
            ) * multiplier
          * varPH_.ddy(wrtVar, cellIndex);
    }
    return returnMe;
}


Foam::tmp<scalarField>
    Foam::admEmpiricalLowerSwitchRateInhibition::ddyFieldFns_VarPH
(
    const admVariable& wrtVar
) const
{
    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            mesh_.nCells(),
            0.0
        )
    );
    scalarField& returnMe(tReturnMe());
    scalarField pH_LL(coeff_pH_LL_.evaluateField());
    scalarField pH_UL(coeff_pH_UL_.evaluateField());
    scalarField pH_var(varPH_.evaluateField());
    scalarField multiplier(mesh_.nCells());
    if (S_Hp_)
    {
        multiplier = -1000.0 * log10(MathConstantScope::e) / pH_var;
        pH_var = - log10(pH_var / 1000.0);
    }
    else
    {
        multiplier = 1.0;
    }
    scalarField denominatorRootRecip
    (
        1.0 / (pH_UL - pH_LL)
    );
    returnMe =
        pos(pH_UL - pH_var) * 6.0 * pH_UL * (pH_var - pH_UL)
      * denominatorRootRecip * denominatorRootRecip
      * exp
        (
            -3.0 * sqr
            (
                (pH_var - pH_UL) * denominatorRootRecip
            )
        ) * multiplier
      * varPH_.ddyField(wrtVar);
    return returnMe;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admEmpiricalLowerSwitchRateInhibition::
    admEmpiricalLowerSwitchRateInhibition
(
    admTime& runTime,
    admVariableManager& admVars,
    admCoefficientManager& admCoeffs,
    const dictionary& dict,
    const word& inhibitionName
)
:
    admRateInhibition(admVars, admCoeffs, inhibitionName),
    
    coeff_pH_LL_
    (
        admCoeffs_.readCoefficient
        (
            dict.subDict(inhibitionName),
            word("pH_LL"),
            word(inhibitionName + "(pH_LL)")
        )
    ),

    coeff_pH_UL_
    (
        admCoeffs_.readCoefficient
        (
            dict.subDict(inhibitionName),
            word("pH_UL"),
            word(inhibitionName + "(pH_UL)")
        )
    ),
    
    // Linked variable can be either pH or S_H+
    varPH_
    (
        admVars_.lookup
        (
            dict.subDict(inhibitionName).found("pH_var")
          ? word(dict.subDict(inhibitionName).lookup("pH_var"))
          : word(dict.subDict(inhibitionName).lookup("S_Hp_var"))
        )
    ),
    
    S_Hp_(dict.subDict(inhibitionName).found("S_Hp_var")),

    ddyFns_(admVars_.nAll()),
    ddyFieldFns_(admVars_.nAll())
{
    if (S_Hp_ && dict.subDict(inhibitionName).found("pH_Var"))
    {
        FatalErrorIn("admEmpiricalLowerSwitchRateInhibition::"
            "admEmpiricalLowerSwitchRateInhibition")
            << "Inhibition " << inhibitionName << "has both pH_var and "
            << "S_Hp_var defined.  Use only the one defined as a standard "
            << "variable in admVariablesDict.  If neither are a standard "
            << "variable, choose only one."
            << abort(FatalError);
    }
    
    if
    (
        (dimensionSet::debug)
     && (
            (coeff_pH_LL_.dimensions() != dimless)
         || (coeff_pH_UL_.dimensions() != dimless)
         || (
                (S_Hp_ && (varPH_.dimensions() != dimMoles / dimVol))
             || (!S_Hp_ && (varPH_.dimensions() != dimMass / dimVol))
            )
        )
    )
    {
        word phVarName("pH_var");
        if (S_Hp_)
        {
            phVarName = "S_Hp_var";
        }
        FatalErrorIn
        (
            "admEmpiricalLowerSwitchRateInhibition"
            "::admEmpiricalLowerSwitchRateInhibition"
        )
            << "Dimension failure for pH-related inhibition " << inhibitionName
            << ".  pH limits and pH_var must be dimensionless, and S_Hp must "
            << "be moles / volume.  Dimensions are:" << endl
            << "pH_LL = " << coeff_pH_LL_.dimensions() << endl
            << "pH_UL = " << coeff_pH_UL_.dimensions() << endl
            << phVarName << " = " << varPH_.dimensions()
            << abort(FatalError);
    }
    
    // Warn if coefficients have non-zero derivatives
    word failedCoeffName(word::null);
    label nzdIndex(0);
    for (; nzdIndex < admVars_.nAll(); nzdIndex++)
    {
        const admVariable& wrtVar(admVars_.all(nzdIndex));
        if (coeff_pH_LL_.ddyNonZero(wrtVar))
        {
            failedCoeffName = coeff_pH_LL_.name();
            break;
        }
        else if (coeff_pH_UL_.ddyNonZero(wrtVar))
        {
            failedCoeffName = coeff_pH_UL_.name();
            break;
        }
    }
    if (failedCoeffName != word::null)
    {
        WarningIn("admEmpiricialLowerSwitchRateInhibition::"
            "admEmpiricialLowerSwitchRateInhibition")
            << "Coefficient " << failedCoeffName << ", specified by "
            << inhibitionName << " rate inhibition, has a non-zero ddy with "
            << "respect to variable " << admVars_.all(nzdIndex).name() << ". "
            << "Calculations of " << inhibitionName << "'s ddy will not take "
            << "this into account.  Use a custom rate inhibition if it is "
            << "required.\n" << endl;
    }
    
    forAll(ddyFns_, index)
    {
        ddyFns_.set
        (
            index,
            new ddyFn
            (
                &Foam::admEmpiricalLowerSwitchRateInhibition::ddyFns_Init
            )
        );
        ddyFieldFns_.set
        (
            index,
            new ddyFieldFn
            (
                &Foam::admEmpiricalLowerSwitchRateInhibition::ddyFieldFns_Init
            )
        );
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::admEmpiricalLowerSwitchRateInhibition::inhibition
(
    label cellIndex
) const
{
    scalar pH_var(varPH_.evaluate(cellIndex));
    if (S_Hp_)
    {
        pH_var = - log10(pH_var / 1000.0);
    }
    scalar pH_LL(coeff_pH_LL_.evaluate(cellIndex));
    scalar pH_UL(coeff_pH_UL_.evaluate(cellIndex));

    return scalar
    (
        pos(pH_UL - pH_var)
      * (
            exp(-3.0 * sqr(pH_var - pH_UL) / (pH_UL - pH_LL)) - 1.0
        ) + 1.0
    );
}


Foam::tmp<scalarField>
    Foam::admEmpiricalLowerSwitchRateInhibition::inhibitionField() const
{
    scalarField pH_var(varPH_.evaluateField());
    if (S_Hp_)
    {
        pH_var = - log10(pH_var / 1000.0);
    }
    scalarField pH_LL(coeff_pH_LL_.evaluateField());
    scalarField pH_UL(coeff_pH_UL_.evaluateField());

    return tmp<scalarField>
    (
        new scalarField
        (
            pos(pH_UL - pH_var)
          * (
                exp(-3 * sqr(pH_var - pH_UL) / (pH_UL - pH_LL)) - 1
            ) + 1
        )
    );
}


bool Foam::admEmpiricalLowerSwitchRateInhibition::inhibitionNonZero() const
{
    return true;
}


Foam::scalar Foam::admEmpiricalLowerSwitchRateInhibition::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    // Call function from jump table
    return (*this.*ddyFns_[wrtVar.globalIndex()])(wrtVar, cellIndex);
}


Foam::tmp<scalarField> Foam::admEmpiricalLowerSwitchRateInhibition::ddyField
(
    const admVariable& wrtVar
) const
{
    // Call function from jump table
    return (*this.*ddyFieldFns_[wrtVar.globalIndex()])(wrtVar);
}


bool Foam::admEmpiricalLowerSwitchRateInhibition::ddyNonZero
(
    const admVariable& wrtVar
) const
{
    return (varPH_.ddyNonZero(wrtVar));
}

// ************************************************************************* //
