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
    Foam::admNonReactingReactionRate

Description
    A reaction rate for the Anaerobic Digester Model (ADM1) that always
    evaluates to zero.  Used to disable reactions.
    
SourceFiles
    admNonReactingReactionRate.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admNonReactingReactionRate_H
#define admNonReactingReactionRate_H

#include "admReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admNonReactingReactionRate Declaration
\*---------------------------------------------------------------------------*/

class admNonReactingReactionRate
:
    public admReactionRate
{
    // Private data

public:

    //- Runtime type information
    TypeName("nonReacting");


    // Constructors

        //- Construct from dictionary, given name
        admNonReactingReactionRate
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            const dictionary& dict,
            const word& reactionName
        );

    // Destructor

        virtual ~admNonReactingReactionRate()
        {}


    // Member functions

    // Member functions

        // Uninhibited value functions
        //  The reactionRate, not including the effects of inhibition

            //- Uninhibited value at a given cell index
            virtual scalar uninhibited(label cellIndex) const
            {
                return scalar(0.0);
            }

            //- Uninhibited value across the entire field
            virtual tmp<scalarField> uninhibitedField() const
            {
                return tmp<scalarField>
                (
                    new scalarField(mesh_.nCells(), scalar(0.0))
                );
            }

        // Uninhibited ddy functions
        //  The reactionRate, not including the effects of inhibition

            //- Uninhibited ddy at a given cell index
            virtual scalar uninhibitedDdy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const
            {
                return scalar(0.0);
            }

            //- Uninhibited ddy across the entire field
            virtual tmp<scalarField> uninhibitedDdyField
            (
                const admVariable& wrtVar
            ) const
            {
                return tmp<scalarField>
                (
                    new scalarField(mesh_.nCells(), scalar(0.0))
                );
            }

            //- True if uninhibited ddy is non-zero
            virtual bool uninhibitedDdyNonZero
            (
                const admVariable& wrtVar
            ) const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
