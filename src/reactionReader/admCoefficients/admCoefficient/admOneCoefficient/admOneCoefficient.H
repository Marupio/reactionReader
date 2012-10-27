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
    Foam::admOneCoefficient

Description
    A coefficient that is one.  Useful for functions returning coefficients,
    when most possibilities are one.
    
SourceFiles
    admOneCoefficient.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admOneCoefficient_H
#define admOneCoefficient_H

#include "admCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admOneCoefficient Declaration
\*---------------------------------------------------------------------------*/

class admOneCoefficient
:
    public admCoefficient
{

public:

    //- Runtime type information
    TypeName("one");


    // Constructors

        //- Build directly (not through runTimeSelection)
        admOneCoefficient
        (
            admVariableManager& admVars,
            const word& coeffName,
            const dimensionSet& dimensions
        );

        //- Construct from dictionary, given name
        admOneCoefficient
        (
            admTime& runTime,
            const fvMesh& mesh,
            admVariableManager& admVars,
            const dictionary& dict,
            const word& coeffName
        );

    // Destructor

        virtual ~admOneCoefficient()
        {}


    // Member functions

        // Evaluation functions

            //- Evaluate a given cell index
            virtual scalar evaluate(label cellIndex) const
            {
                return 1.0;
            }

            //- Evaluate the entire internal field
            virtual tmp<scalarField> evaluateField() const
            {
                return tmp<scalarField>
                (
                    new scalarField(admVars_.mesh().nCells(), 1.0)
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
