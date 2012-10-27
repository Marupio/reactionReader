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
    Foam::admNonCompetitiveRateInhibition

Description
    A reaction rate inhibition multiplier for the Anaerobic Digester Model
    (ADM1) given by:

    \f[
        I = \frac{K_I}{K_I + S_I}
    \f]

    where:
        \f$ K_I \f$ is a numerical constant; and
        \f$ S_I \f$ is the variable it depends on.
    
SourceFiles
    admNonCompetitiveRateInhibition.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admNonCompetitiveRateInhibition_H
#define admNonCompetitiveRateInhibition_H

#include "admRateInhibition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class admNonCompetitiveRateInhibition Declaration
\*---------------------------------------------------------------------------*/

class admNonCompetitiveRateInhibition
:
    public admRateInhibition
{
    // Private data

        // Function pointers
        //  For efficiency we use a jump table for ddy and ddyField functions

            //- ddy function typedef
            typedef scalar (Foam::admNonCompetitiveRateInhibition::*ddyFn)
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- ddyField function typedef
            typedef tmp<scalarField>
            (Foam::admNonCompetitiveRateInhibition::*ddyFieldFn)
            (
                const admVariable& wrtVar
            ) const;


        //- Associated coefficients
        declareAdmJumpTableCoefficient(coeff_KI_);

        //- Associated variables
        declareAdmJumpTableVariable(var_SI_);

    // Private member functions
    
        // Universal function pointers
        
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

public:

    //- Runtime type information
    TypeName("nonCompetitive");


    // Constructors

        //- Construct from dictionary, given name
        admNonCompetitiveRateInhibition
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            const dictionary& dict,
            const word& inhibitionName
        );

    // Destructor

        virtual ~admNonCompetitiveRateInhibition()
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
