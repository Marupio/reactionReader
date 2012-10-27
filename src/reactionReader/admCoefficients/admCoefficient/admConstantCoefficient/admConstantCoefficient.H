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
    Foam::admConstantCoefficient

Description
    A time and space invariant (fully constant) coefficient for the Anaerobic
    Digester Model (ADM1)
    
SourceFiles
    admConstantCoefficient.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admConstantCoefficient_H
#define admConstantCoefficient_H

#include "admCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class admConstantCoefficient Declaration
\*---------------------------------------------------------------------------*/

class admConstantCoefficient
:
    public admCoefficient
{

private:

    //- Coefficient value (dimensionedScalar components are held by base class)
    mutable scalar value_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary, given name
        admConstantCoefficient
        (
            admTime& runTime,
            const fvMesh& mesh,
            admVariableManager& admVars,
            const dictionary& dict,
            const word& coeffName
        );

    // Destructor

        virtual ~admConstantCoefficient()
        {}


    // Member functions

        // Evaluation functions

            virtual scalar evaluate(label cellIndex) const
            {
                return value_;
            }

            virtual tmp<scalarField> evaluateField() const
            {
                return tmp<scalarField>
                (
                    new scalarField
                    (
                        admVars_.mesh().nCells(),
                        value_
                    )
                );
            }

            //- Returns true if the coefficient does not change spatially
            virtual bool uniform() const
            {
                return true;
            }

            //- Returns true if the coefficient is invariant with all variables
            virtual bool constant() const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
