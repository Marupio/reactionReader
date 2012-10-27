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
    Foam::admSimpleGasReactionRate

Description
    A reaction rate for the Anaerobic Digester Model (ADM1) given by:

    \f[
        \rho = k_La \left ( S_{var} + n \ K_H p_{var} \right )
    \f]

    where:
        \f$ k_La \f$ is the overall gas transfer coefficient;
        \f$ S_{var} \f$ is the variable for the dissolved gas component;
        n is the fraction of COD equivalents of the gas;
        \f$ K_H \f$ is the coefficient of Henry's Law; and
        \f$ p_{var} \f$ is the variable for the gas partial pressure.
        
SourceFiles
    admSimpleGasReactionRate.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admSimpleGasReactionRate_H
#define admSimpleGasReactionRate_H

#include "admReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admSimpleGasReactionRate Declaration
\*---------------------------------------------------------------------------*/

class admSimpleGasReactionRate
:
    public admReactionRate
{
    // Private data

        // Function pointers
        //  For efficiency we use a set of jump tables for the ddy functions

            //- ddy function typedef
            typedef scalar (Foam::admSimpleGasReactionRate::*ddyFn)
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- ddyField function typedef
            typedef tmp<scalarField>
            (Foam::admSimpleGasReactionRate::*ddyFieldFn)
            (
                const admVariable& wrtVar
            ) const;

        //- Associated coefficients
        declareAdmJumpTableCoefficient(coeff_kLa_);
        declareAdmJumpTableCoefficient(coeff_n_);
        declareAdmJumpTableCoefficient(coeff_KH_);

        //- Associated variables
        declareAdmJumpTableVariable(var_S_);
        declareAdmJumpTableVariable(var_P_);

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
    TypeName("simpleGas");

    // Constructors

        //- Construct from dictionary, given name
        admSimpleGasReactionRate
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            const dictionary& dict,
            const word& reactionName
        );

    // Destructor

        virtual ~admSimpleGasReactionRate()
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
