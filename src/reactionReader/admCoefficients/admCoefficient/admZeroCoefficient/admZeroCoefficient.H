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
    Foam::admZeroCoefficient

Description
    A coefficient that is zero.  Useful for functions returning coefficients,
    when most possibilities are zero.
    
SourceFiles
    admZeroCoefficient.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admZeroCoefficient_H
#define admZeroCoefficient_H

#include "admCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admZeroCoefficient Declaration
\*---------------------------------------------------------------------------*/

class admZeroCoefficient
:
    public admCoefficient
{

public:

    //- Runtime type information
    TypeName("zero");


    // Constructors

        //- Build directly (not through runTimeSelection)
        admZeroCoefficient
        (
            admVariableManager& admVars,
            const word& coeffName,
            const dimensionSet& dimensions
        );

        //- Construct from dictionary, given name
        admZeroCoefficient
        (
            admTime& runTime,
            const fvMesh& mesh,
            admVariableManager& admVars,
            const dictionary& dict,
            const word& coeffName
        );

    // Destructor

        virtual ~admZeroCoefficient()
        {}


    // Member functions

        // Evaluation functions

            //- Evaluate a given cell index
            virtual scalar evaluate(label cellIndex) const
            {
                return 0.0;
            }

            //- Evaluate the entire internal field
            virtual tmp<scalarField> evaluateField() const
            {
                return tmp<scalarField>
                (
                    new scalarField(admVars_.mesh().nCells(), 0.0)
                );
            }

            //- True if evaluate is non-zero
            virtual bool evaluateNonZero() const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
