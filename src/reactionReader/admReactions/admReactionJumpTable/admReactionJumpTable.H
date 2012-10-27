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

Class
    Foam::admReactionJumpTable

Description
    A set of macros to facilitate the use of jump tables in the derived classes
    for admReactionReactionRate and admRateInhibition.  Tests show a 10%
    improvement in overall model simulation time.
    
SourceFiles

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/


#ifndef admReactionJumpTable_H
#define admReactionJumpTable_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// private macro

#define declareAdmReactionJumpTables(coeffOrVarName)                          \
    /* Jump tables */                                                         \
    mutable PtrList<ddyFn> ddyTable##coeffOrVarName;                          \
    mutable PtrList<ddyFieldFn> ddyFieldTable##coeffOrVarName;                \
                                                                              \
    /* Initialization functions */                                            \
    scalar ddyFns_init##coeffOrVarName                                        \
    (                                                                         \
        const admVariable& wrtVar,                                            \
        label cellIndex                                                       \
    ) const;                                                                  \
    tmp<scalarField> ddyFieldFns_init##coeffOrVarName                         \
    (                                                                         \
        const admVariable& wrtVar                                             \
    ) const;                                                                  \
                                                                              \
    /* Partial derivative functions */                                        \
    scalar ddyFns##coeffOrVarName                                             \
    (                                                                         \
        const admVariable& wrtVar,                                            \
        label cellIndex                                                       \
    ) const;                                                                  \
    tmp<scalarField> ddyFieldFns##coeffOrVarName                              \
    (                                                                         \
        const admVariable& wrtVar                                             \
    ) const

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in header file
// for adding a jump-table-enabled variable

#define declareAdmJumpTableVariable(variableName)                             \
    /* Variable reference */                                                  \
    const admVariable& variableName;                                          \
    declareAdmReactionJumpTables(variableName)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in header file
// for adding a jump-table-enabled coefficient

#define declareAdmJumpTableCoefficient(coefficientName)                       \
    /* Variable reference */                                                  \
    const admCoefficient& coefficientName;                                    \
    declareAdmReactionJumpTables(coefficientName)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file
// for adding jump table initialization functions

#define defineAdmJumpTableCoefficientOrVariable(coeffOrVarName, className)    \
Foam::scalar Foam::className::ddyFns_init##coeffOrVarName                     \
(                                                                             \
    const admVariable& wrtVar,                                                \
    label cellIndex                                                           \
) const                                                                       \
{                                                                             \
    label globalIndex(wrtVar.globalIndex());                                  \
    if (coeffOrVarName.ddyNonZero(wrtVar))                                    \
    {                                                                         \
        ddyTable##coeffOrVarName[globalIndex] =                               \
            &Foam::className::ddyFns##coeffOrVarName;                         \
        ddyFieldTable##coeffOrVarName[globalIndex] =                          \
            &Foam::className::ddyFieldFns##coeffOrVarName;                    \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        ddyTable##coeffOrVarName[globalIndex] =                               \
            &Foam::className::ddyFns_0;                                       \
        ddyFieldTable##coeffOrVarName[globalIndex] =                          \
            &Foam::className::ddyFieldFns_0;                                  \
    }                                                                         \
    /* Call the originally requested function */                              \
    return (*this.*ddyTable##coeffOrVarName[globalIndex])(wrtVar, cellIndex); \
}                                                                             \
                                                                              \
                                                                              \
Foam::tmp<scalarField> Foam::className::ddyFieldFns_init##coeffOrVarName      \
(                                                                             \
    const admVariable& wrtVar                                                 \
) const                                                                       \
{                                                                             \
    ddyFns_init##coeffOrVarName(wrtVar, 0);                                   \
    /* Call the originally requested function */                              \
    return                                                                    \
        (*this.*ddyFieldTable##coeffOrVarName[wrtVar.globalIndex()])(wrtVar); \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file
// for adding point derivative functions with respect to the given coefficient
// or variable.  Replaces the return type, function name, and parameter list.
// User must define the function implementation beneath, in braces {}.

#define admJumpTableDerivative(coeffOrVarName, className)                     \
Foam::scalar Foam::className::ddyFns##coeffOrVarName                          \
(                                                                             \
    const admVariable& wrtVar,                                                \
    label cellIndex                                                           \
) const

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file
// for adding field derivative functions with respect to the given coefficient
// or variable.  Replaces the return type, function name, and parameter list.
// User must define the function implementation beneath, in braces {}.

#define admJumpTableFieldDerivative(coeffOrVarName, className)                \
Foam::tmp<scalarField> Foam::className::ddyFieldFns##coeffOrVarName           \
(                                                                             \
    const admVariable& wrtVar                                                 \
) const

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file
// for calling a point derivative function from the jump table

#define admJumpTableRetrieveDerivative(coeffOrVarName)                        \
    (*this.*ddyTable##coeffOrVarName[wrtVar.globalIndex()])(wrtVar, cellIndex)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file
// for calling a field derivative function from the jump table

#define admJumpTableRetrieveFieldDerivative(coeffOrVarName)                   \
    (*this.*ddyFieldTable##coeffOrVarName[wrtVar.globalIndex()])(wrtVar)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file for admReactionRates
// use in constructor initialization table to initialize a coefficient for a
// derived admReactionRate

#define initializeAdmReactionRateCoefficient(coefficientName, lookupWord)     \
    coefficientName                                                           \
    (                                                                         \
        admCoeffs_.readCoefficient                                            \
        (                                                                     \
            dict.subDict(reactionName).subDict("rate"),                       \
            word(#lookupWord),                                                \
            word(reactionName + "(rate(" #lookupWord "))")                    \
        )                                                                     \
    ),                                                                        \
    ddyTable##coefficientName(admVars_.nAll()),                               \
    ddyFieldTable##coefficientName(admVars_.nAll())

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file for admRateInhibitions
// use in constructor initialization table to initialize a coefficient for a
// derived admRateInhibition

#define initializeAdmRateInhibitionCoefficient(coefficientName, lookupWord)   \
    coefficientName                                                           \
    (                                                                         \
        admCoeffs_.readCoefficient                                            \
        (                                                                     \
            dict.subDict(inhibitionName),                                     \
            word(#lookupWord),                                                \
            word(inhibitionName + "(" #lookupWord ")")                        \
        )                                                                     \
    ),                                                                        \
    ddyTable##coefficientName(admVars_.nAll()),                               \
    ddyFieldTable##coefficientName(admVars_.nAll())

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file for admReactionRates
// use in constructor initialization table to initialize a variable for a
// derived admReactionRate

#define initializeAdmReactionRateVariable(variableName, lookupWord)           \
    variableName                                                              \
    (                                                                         \
        admVars_.lookup                                                       \
        (                                                                     \
            word                                                              \
            (                                                                 \
                dict.subDict(reactionName).subDict("rate")                    \
                    .lookup(#lookupWord)                                      \
            )                                                                 \
        )                                                                     \
    ),                                                                        \
    ddyTable##variableName(admVars_.nAll()),                                  \
    ddyFieldTable##variableName(admVars_.nAll())

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file for admRateInhibitions
// use in constructor initialization table to initialize a variable for a
// derived admRateInhibition

#define initializeAdmRateInhibitionVariable(variableName, lookupWord)         \
    variableName                                                              \
    (                                                                         \
        admVars_.lookup                                                       \
        (                                                                     \
            word                                                              \
            (                                                                 \
                dict.subDict(inhibitionName).lookup(#lookupWord)              \
            )                                                                 \
        )                                                                     \
    ),                                                                        \
    ddyTable##variableName(admVars_.nAll()),                                  \
    ddyFieldTable##variableName(admVars_.nAll())

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// public macro - used in source file
// use in constructor body to allocate memory for jump table components

#define allocateAdmCoefficientOrVariable(coeffOrVarName, className)           \
    for (label i(0); i < admVars_.nAll(); i++)                                \
    {                                                                         \
        ddyTable##coeffOrVarName.set                                          \
        (                                                                     \
            i,                                                                \
            new ddyFn(&Foam::className::ddyFns_init##coeffOrVarName)          \
        );                                                                    \
        ddyFieldTable##coeffOrVarName.set                                     \
        (                                                                     \
            i,                                                                \
            new ddyFieldFn                                                    \
            (                                                                 \
                &Foam::className::ddyFieldFns_init##coeffOrVarName            \
            )                                                                 \
        );                                                                    \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
