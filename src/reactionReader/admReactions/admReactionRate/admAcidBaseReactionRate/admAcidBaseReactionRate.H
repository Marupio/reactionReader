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
    Foam::admAcidBaseReactionRate

Description
    A reaction rate for the Anaerobic Digester Model (ADM1) given by:

    \f[
        \rho = k_{AB} \left ( S_{var^-} S_{H^+} - K_a S_{hvar} \right )
    \f]
    
    or

    \f[
        \rho = k_{AB} 
            \left \lbrack 
                S_{var^-} 
                \left (
                    K_a + S_{H^+}
                \right )
                - K_a S_{var}
            \right \rbrack
    \f]
    

    where:
        \f$ k_{AB} \f$ is the acid-base kinetic parameter;
        \f$ S_{var^-} \f$ is the concentration of the compound's negative ion
            form;
        \f$ S_{H^+} \f$ is the concentration of hydrogen ions;
        \f$ K_{a} \f$ is the acid-base equilibrium coefficient;
        \f$ S_{hvar} \f$ is the concentration of the compound's acid form; and
        \f$ S_{hvar} \f$ is the total concentration of the compound.
        
SourceFiles
    admAcidBaseReactionRate.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admAcidBaseReactionRate_H
#define admAcidBaseReactionRate_H

#include "admReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admAcidBaseReactionRate Declaration
\*---------------------------------------------------------------------------*/

class admAcidBaseReactionRate
:
    public admReactionRate
{
    // Private data

        // Function pointers
        //  For efficiency we use a jump table for ddy and ddyField functions

            //- ddy function typedef
            typedef scalar (Foam::admAcidBaseReactionRate::*ddyFn)
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- ddyField function typedef
            typedef tmp<scalarField>
            (Foam::admAcidBaseReactionRate::*ddyFieldFn)
            (
                const admVariable& wrtVar
            ) const;

        //- Associated coefficients
        declareAdmJumpTableCoefficient(coeff_kAB_);
        declareAdmJumpTableCoefficient(coeff_Ka_);

        //- Associated variables
        declareAdmJumpTableVariable(var_Sm_);
        declareAdmJumpTableVariable(var_S_);
        declareAdmJumpTableVariable(var_Shp_);

        //- 1.0 if S_var is actually S_hvar, 0.0 otherwise
        scalar using_S_hvar_;

    // Private member functions
    
        // Universal function pointers
        
            //- Returns zero
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
    TypeName("acidBase");

    // Constructors

        //- Construct from dictionary, given name
        admAcidBaseReactionRate
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            const dictionary& dict,
            const word& reactionName
        );

    // Destructor

        virtual ~admAcidBaseReactionRate()
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
