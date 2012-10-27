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
    Foam::admCompetitiveUptakeRateInhibition

Description
    A reaction rate inhibition multiplier for the Anaerobic Digester Model
    (ADM1) given by:

    \f[
        I = \frac{S_{var}}{S_{var} + S_{I}}
    \f]

    where:
        \f$ S_{var} \f$ is the substrate variable for the inhibited process;
            and
        \f$ S_I \f$ is the inhibitor variable.
    
SourceFiles
    admCompetitiveUptakeRateInhibition.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admCompetitiveUptakeRateInhibition_H
#define admCompetitiveUptakeRateInhibition_H

#include "admRateInhibition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class admCompetitiveUptakeRateInhibition Declaration
\*---------------------------------------------------------------------------*/

class admCompetitiveUptakeRateInhibition
:
    public admRateInhibition
{
    // Private data
        
        // Function pointers
        //  For efficiency we use a jump table for ddy and ddyField functions

            //- ddy function typedef
            typedef scalar (Foam::admCompetitiveUptakeRateInhibition::*ddyFn)
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- ddyField function typedef
            typedef tmp<scalarField>
            (Foam::admCompetitiveUptakeRateInhibition::*ddyFieldFn)
            (
                const admVariable& wrtVar
            ) const;

        //- Associated variables
        
            //- S_var - the process substrate
            declareAdmJumpTableVariable(var_S_);

            //- S_I - the inhibitor
            declareAdmJumpTableVariable(var_SI_);
            
        //- Other user-defined settings
        
            //- Stabilise offset, a small quantity added to the denominator to
            //  prevent division by zero
            scalar stabiliseEpsilon_;

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
    TypeName("competitiveUptake");

    // Constructors

        //- Construct from dictionary, given name
        admCompetitiveUptakeRateInhibition
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            const dictionary& dict,
            const word& inhibitionName
        );

    // Destructor

        virtual ~admCompetitiveUptakeRateInhibition()
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
