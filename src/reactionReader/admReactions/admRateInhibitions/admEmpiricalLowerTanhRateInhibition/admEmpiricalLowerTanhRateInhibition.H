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
    Foam::admEmpiricalLowerTanhRateInhibition

Description
    A reaction rate inhibition multiplier for the Anaerobic Digester Model
    (ADM1) given by:

    \f[
    I =
    \frac{1}{2}
    \left \lbrack
        1 + \tanh
        \left (
            a \frac{pH_{var}}{pH_{avg}} - 1
        \right )
    \right \rbrack
    \f]

    where:
        a is a numerical constant;
        \f$ pH_{avg} \f$ is given by:
    
    \f[
    pH_{avg} = \frac{pH_{LL} + pH_{UL}}{2}
    \f]
    
    and:
        \f$ pH_{LL} \f$ is the lower pH limit;
        \f$ pH_{UL} \f$ is the upper pH limit; and
        \f$ pH_{var} \f$ is the variable pH.
        
SourceFiles
    admEmpiricalLowerTanhRateInhibition.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admEmpiricalLowerTanhRateInhibition_H
#define admEmpiricalLowerTanhRateInhibition_H

#include "admRateInhibition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class admEmpiricalLowerTanhRateInhibition Declaration
\*---------------------------------------------------------------------------*/

class admEmpiricalLowerTanhRateInhibition
:
    public admRateInhibition
{
    // Private data

        //- Numerical constants (coefficient names)
        const admCoefficient& coeff_a_;
        const admCoefficient& coeff_pH_LL_;
        const admCoefficient& coeff_pH_UL_;
        
        //- Associated variable - pH
        const admVariable& varPH_;

        //- true if the associated variable is S_H+ instead of pH
        bool S_Hp_;

        // Function pointers
        //  For efficiency we use a jump table for ddy and ddyField functions

            //- ddy function typedef
            typedef scalar (Foam::admEmpiricalLowerTanhRateInhibition::*ddyFn)
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- ddyField function typedef
            typedef tmp<scalarField>
            (Foam::admEmpiricalLowerTanhRateInhibition::*ddyFieldFn)
            (
                const admVariable& wrtVar
            ) const;
            
            //- ddy jump table
            mutable PtrList<ddyFn> ddyFns_;

            //- ddyField jump table
            mutable PtrList<ddyFieldFn> ddyFieldFns_;

    // Private member functions
    
        // Function pointers
        
            //- Initialize function pointers
            scalar ddyFns_Init
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            tmp<scalarField> ddyFieldFns_Init
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

            // Used when wrt variable ddy is nonzero to varPH_
            scalar ddyFns_VarPH
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            tmp<scalarField> ddyFieldFns_VarPH
            (
                const admVariable& wrtVar
            ) const;

public:

    //- Runtime type information
    TypeName("empiricalLowerTanh");

    // Constructors

        //- Construct from dictionary, given name
        admEmpiricalLowerTanhRateInhibition
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            const dictionary& dict,
            const word& inhibitionName
        );

    // Destructor

        virtual ~admEmpiricalLowerTanhRateInhibition()
        {}

    // Member functions
    
        // Inhibition functions

            //- Evaluate inhibition at a given cell index
            virtual scalar inhibition(label cellIndex) const;

            //- Evaluate inhibition throughout entire internal field
            virtual tmp<scalarField> inhibitionField() const;

            //- True if inhibition is non-zero
            virtual bool inhibitionNonZero() const;

        // Ddy functions

            virtual scalar ddy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            virtual tmp<scalarField> ddyField
            (
                const admVariable& wrtVar
            ) const;

            virtual bool ddyNonZero(const admVariable& wrtVar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
