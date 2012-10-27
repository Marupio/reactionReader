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
    Foam::admStandardVariable

Description
    Anaerobic digester standard variable - this is the "main" variable type in
    the model, and is intended to be solved using conventional methods.

SourceFiles
    admStandardVariableI.H
    admStandardVariable.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admStandardVariable_H
#define admStandardVariable_H

#include "admVariable.H"
#include "equation.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class admStandardVariable Declaration
\*---------------------------------------------------------------------------*/

class admStandardVariable
:
    public admVariable
{

    // Private member data

        //- Pointer to the volScalarField containing the data
        autoPtr<volScalarField> dataPtr_;

        //- Diffusion coefficient
        dimensionedScalar gamma_;

        //- Main loop convergence
        scalar convergence_;
        
        //- Minimum error scaling factor
        scalar minErrorScale_;
        
        //- When true, a user-defined function changes this variable and must
        //  be accounted for in the main solver
        Switch changedByUdf_;

private:

    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        admStandardVariable(const admStandardVariable&);
        void operator=(const admStandardVariable&);

public:

    TypeName("admStandardVariable");

    // Constructors

        //- Construct null - internal use only - to permit returning null
        admStandardVariable
        (
            admTime& runTime,
            admVariableManager& admVars
        );

        //- Construct from dictionary
        admStandardVariable
        (
            admTime& runTime,
            admVariableManager& admVars,
            const dictionary& dict,
            const label globalIndex,
            const label localIndex,
            const word& name
        );

    // Destructor

        virtual ~admStandardVariable()
        {}


    // Member Functions

        // Access
        
            //- Access the diffusion coefficient
            inline const dimensionedScalar& gamma() const;
            inline dimensionedScalar& gamma();
            
            //- Return main loop convergence
            inline const scalar& convergence() const;
            
            //- Return minErrorScale
            inline const scalar& minErrorScale() const;
        
            //- Return changedByUdf
            inline const Switch& changedByUdf() const;

            //- Access the data directly
            inline const volScalarField& operator()() const;
            inline volScalarField& operator()();

            //- Access the data directly by cellIndex
            inline const scalar& operator()(label cellIndex) const;
            inline scalar& operator()(label cellIndex);

        // Access - semantic shortcuts
        
            //- Variable name
            virtual const word& name() const;

            //- Variable dimensions
            virtual const dimensionSet& dimensions() const;

        // Assign functions

            //- Assign a value at a given cell index
            virtual void assign(const scalar& scalarIn, label cellIndex);

            //- Assign the dimensions - only checks dimensions
            virtual void assignDims
            (
                const dimensionSet& dimsIn
            );
            
            //- Assign a dimensioned value at a given cell index
            virtual void assignDimensioned
            (
                const dimensionedScalar& dScalarIn,
                label cellIndex
            );

            //- Assign values to the entire field
            virtual void assignField
            (
                const scalarField& sfIn
            );

            //- Assign dimensioned values to the entire field            
            virtual void assignDimensionedField
            (
                const dimensionedScalarField& dsfIn
            );

        // Evaluate functions
        //  Evaluate functions for variables return const references where
        //  possible to prevent redundant copying, as opposed to straight up
        //  objects passed in the rest of the framework

            //- Evaluate at a given cell index
            virtual const scalar& evaluate(label cellIndex) const;

            //- Return the name associated with an evaluation
            virtual word evaluateName() const;

            //- Return the dimensions after evaluation
            virtual const dimensionSet& evaluateDims() const;

            //- Evaluate at a given cell index with dimensions
            virtual dimensionedScalar evaluateDimensioned
            (
                label cellIndex
            ) const;

            //- Evaluate the entire internal field
            virtual const scalarField& evaluateField() const;

            //- Evaluate the entire internal field with dimensions
            virtual const dimensionedScalarField&
                evaluateDimensionedField() const;

            //- Evaluate the entire Geometric field
            virtual const volScalarField& evaluateGeometricField() const;

            //- True if evaluate is non-zero
            virtual bool evaluateNonZero() const;

        // Ddy functions

            //- Derivative with respect to a given variable at a given cell
            //  index
            virtual scalar ddy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- Return the name associated with the derivative with respect to
            //  a given variable
            virtual word ddyName
            (
                const admVariable& wrtVar
            ) const;

            //- Return the dimensions the derivative with respect to a given
            //  variable
            virtual dimensionSet ddyDims
            (
                const admVariable& wrtVar
            ) const;

            //- Derivative with respect to a given variable at a given cell
            //  index with dimensions
            virtual dimensionedScalar ddyDimensioned
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- Derivative with respect to a given variable across the entire
            //  internal field
            virtual tmp<scalarField> ddyField
            (
                const admVariable& wrtVar
            ) const;

            //- Derivative with respect to a given variable across the entire
            //  internal field with dimensions
            virtual tmp<dimensionedScalarField> ddyDimensionedField
            (
                const admVariable& wrtVar
            ) const;

            //- True if derivative with respect to a given variable is non-zero
            virtual bool ddyNonZero(const admVariable& wrtVar) const;

        // Housekeeping
            
            //- Correct boundary conditions
            virtual void correctBoundaryConditions();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admStandardVariableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
