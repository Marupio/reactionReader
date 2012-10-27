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
    Foam::admCompetitiveMonodReactionRate

Description
    A reaction rate for the Anaerobic Digester Model (ADM1) given by:

    \f[
        \rho = \frac
        {
            k_m K_I S_{var} X_{var}
        }
        {
            K_s \left ( K_I + S_I \right ) + S_{var} K_I
        }
    \f]

    where:
    
    \f$ k_m \f$ is the kinetic rate parameter;
    \f$ K_s \f$ is the half-saturation value of the compound;
    \f$ K_I \f$ is the half-saturation value of the inhibitor;
    \f$ S_{var} \f$ is the soluble dissolved concentration of the compound;
    \f$ S_I \f$ is the soluble dissolved concentration of the inhibitor; and
    \f$ X_{var} \f$ is the particulate concentration of the compound.

SourceFiles
    admCompetitiveMonodReactionRate.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admCompetitiveMonodReactionRate_H
#define admCompetitiveMonodReactionRate_H

#include "admReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admCompetitiveMonodReactionRate Declaration
\*---------------------------------------------------------------------------*/

class admCompetitiveMonodReactionRate
:
    public admReactionRate
{
    // Private data

        // Function pointers
        //  For efficiency we use a jump table for ddy and ddyField functions

            //- ddy function typedef
            typedef scalar (Foam::admCompetitiveMonodReactionRate::*ddyFn)
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- ddyField function typedef
            typedef tmp<scalarField>
            (Foam::admCompetitiveMonodReactionRate::*ddyFieldFn)
            (
                const admVariable& wrtVar
            ) const;
            
        //- Associated coefficients
        declareAdmJumpTableCoefficient(coeff_km_);
        declareAdmJumpTableCoefficient(coeff_Ks_);
        declareAdmJumpTableCoefficient(coeff_KI_);

        //- Associated variables
        declareAdmJumpTableVariable(var_S_);
        declareAdmJumpTableVariable(var_X_);
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
    TypeName("competitiveMonod");


    // Constructors

        //- Construct from dictionary, given name
        admCompetitiveMonodReactionRate
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            const dictionary& dict,
            const word& reactionName
        );

    // Destructor

        virtual ~admCompetitiveMonodReactionRate()
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
