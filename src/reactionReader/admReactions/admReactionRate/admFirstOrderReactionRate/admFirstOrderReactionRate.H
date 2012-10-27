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
    Foam::admFirstOrderReactionRate

Description
    A reaction rate for the Anaerobic Digester Model (ADM1) given by:

    \f[
        \rho = k * var
    \f]

    where:
        k is the first order numerical constant; and
        var is the variable it depends on.
    
SourceFiles
    admFirstOrderReactionRate.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admFirstOrderReactionRate_H
#define admFirstOrderReactionRate_H

#include "admReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admFirstOrderReactionRate Declaration
\*---------------------------------------------------------------------------*/

class admFirstOrderReactionRate
:
    public admReactionRate
{
    // Private data

        //- Associated coefficient
        const admCoefficient& coeff_k_;

        //- Associated variable
        const admVariable& var_;

        // Function pointers
        //  For efficiency we use a set of jump tables for the ddy functions

            //- ddy function typedef
            typedef scalar (Foam::admFirstOrderReactionRate::*ddyFn)
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- ddyField function typedef
            typedef tmp<scalarField>
            (Foam::admFirstOrderReactionRate::*ddyFieldFn)
            (
                const admVariable& wrtVar
            ) const;
            
            //- ddy jump tables
            mutable PtrList<ddyFn> ddyTable_k_;
            mutable PtrList<ddyFn> ddyTable_var_;

            //- ddyField jump tables
            mutable PtrList<ddyFieldFn> ddyFieldTable_k_;
            mutable PtrList<ddyFieldFn> ddyFieldTable_var_;

    // Private member functions
    
        // Function pointers
        
            //- Initialize function pointers
            scalar ddyFns_init_k
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            scalar ddyFns_init_var
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            tmp<scalarField> ddyFieldFns_init_k
            (
                const admVariable& wrtVar
            ) const;
            tmp<scalarField> ddyFieldFns_init_var
            (
                const admVariable& wrtVar
            ) const;
            
            //- Used when wrt variable ddy is zero
            scalar ddyFns_0
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            tmp<scalarField> ddyFieldFns_0
            (
                const admVariable& wrtVar
            ) const;

            //- Used when wrt variable ddy is nonzero to k_
            scalar ddyFns_k
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            tmp<scalarField> ddyFieldFns_k
            (
                const admVariable& wrtVar
            ) const;

            //- Used when wrt variable ddy is nonzero to var_
            scalar ddyFns_var
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            tmp<scalarField> ddyFieldFns_var
            (
                const admVariable& wrtVar
            ) const;

public:

    //- Runtime type information
    TypeName("firstOrder");


    // Constructors

        //- Construct from dictionary, given name
        admFirstOrderReactionRate
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            const dictionary& dict,
            const word& reactionName
        );

    // Destructor

        virtual ~admFirstOrderReactionRate()
        {}

    // Member functions

        // Uninhibited value functions
        //  The reactionRate, not including the effects of inhibition

            //- Uninhibited value at a given cell index
            virtual scalar uninhibited(label cellIndex) const;

            //- Uninhibited value across the entire field
            virtual tmp<scalarField> uninhibitedField() const;

        // Uninhibited ddy functions
        //  The reactionRate, not including the effects of inhibition

            //- Uninhibited ddy at a given cell index
            virtual scalar uninhibitedDdy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- Uninhibited ddy across the entire field
            virtual tmp<scalarField> uninhibitedDdyField
            (
                const admVariable& wrtVar
            ) const;

            //- True if uninhibited ddy is non-zero
            virtual bool uninhibitedDdyNonZero
            (
                const admVariable& wrtVar
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
