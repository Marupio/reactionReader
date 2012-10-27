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
    Foam::admCoefficient

Description
    Abstract base class for coefficients in the anaerobic digester model, ADM1.
    This class contains the interface for all coefficients.
    
SourceFiles
    admCoefficient.C
    newadmCoefficient.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admCoefficient_H
#define admCoefficient_H

#include "refCount.H"
#include "autoPtr.H"
#include "dimensionedScalar.H"
#include "admTime.H"
#include "fvMesh.H"
#include "admVariableManager.H"
#include "admCalculusInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class admCoefficient Declaration
\*---------------------------------------------------------------------------*/

class admCoefficient
:
    public admCalculusInterface
{

    // Private member functions
    
        //- Disallow default bitwise copy construct and assignment
        void operator=(const admCoefficient&);
        admCoefficient(const admCoefficient&);

public:

    //- Runtime type information
    TypeName("admCoefficient");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            admCoefficient,
            dictionary,
            (
                admTime& runTime,
                const fvMesh& mesh,
                admVariableManager& admVars,
                const dictionary& dict,
                const word& coeffName
            ),
            (runTime, mesh, admVars, dict, coeffName)
        );


    // Constructors

        //- Construct from components
        admCoefficient
        (
            admVariableManager& admVars,
            const word& name,
            const dimensionSet& dimensions
        );

    // Selectors


        //- Return a reference to a new admCoefficient
        //  *** Located in newAdmCoefficient.C ***
        static autoPtr<admCoefficient> New
        (
            admTime& runTime,
            const fvMesh& mesh,
            admVariableManager& admVars,
            const dictionary& dict,
            const word& coeffName
        );


    // Destructor

        virtual ~admCoefficient()
        {}


    // Member functions

        // Access

            //- Returns true if the coefficient does not change spatially
            virtual bool uniform() const
            {
                return false;
            }
            
            //- Returns true if the coefficient is invariant with all variables
            virtual bool constant() const
            {
                return false;
            }


        // ddt functions

            //- Evaluate ddt at a given cell index
            virtual scalar ddt(label cellIndex) const
            {
                return 0;
            }

            //- Return the name associated with a ddt evaluation
            word ddtName() const
            {
                return word("ddt(" + name() + ")");
            }

            //- Return the dimensions after a ddt evaluation
            dimensionSet ddtDims() const
            {
                return dimensions() / dimTime;
            }

            //- Evaluate ddt at a given cell index with dimensions
            dimensionedScalar ddtDimensioned
            (
                label cellIndex
            ) const;

            //- Evaluate ddt over the entire internal field
            virtual tmp<scalarField> ddtField() const;

            //- Evaluate ddt over the entire internal field with dimensions
            tmp<dimensionedScalarField> ddtDimensionedField();

            //- True if ddt is non-zero
            virtual bool ddtNonZero() const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
