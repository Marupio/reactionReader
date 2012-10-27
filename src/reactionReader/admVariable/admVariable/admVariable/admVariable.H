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
    Foam::admVariable

Description
    Anaerobic digester variable, base class.

SourceFiles
    admVariableI.H
    admVariable.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admVariable_H
#define admVariable_H

#include "admTime.H"
#include "word.H"
#include "dimensionedScalar.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class admVariableManager;

/*---------------------------------------------------------------------------*\
                           Class admVariable Declaration
\*---------------------------------------------------------------------------*/

class admVariable
:
    public dictionary
{

public:
    enum varTypeEnum
    {
        vtnone,
        vtstandard,
        vtimplicit,
        vtderived
    };

protected:

    //- Reference to admTime
    admTime& runTime_;

    //- Reference to the mesh
    const fvMesh& mesh_;
    
    //- Reference to variable manager
    admVariableManager& admVars_;

    //- Base index - position in admVariableManager::all_;
    const label globalIndex_;

    //- Type index - position in admVariableManager::[derivedType]_
    //  e.g. standard_, implicit_, or derived_;
    const label localIndex_;

    //- Type of variable - set by derived classes
    varTypeEnum varType_;
    
private:

    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        admVariable(const admVariable&);
        void operator=(const admVariable&);

public:

    TypeName("admVariable");

    // Constructors

        //- Construct from dictionary
        admVariable
        (
            admTime& runTime,
            admVariableManager& admVars,
            const dictionary& dict,
            const label globalIndex,
            const label localIndex
        );

    // Destructor

        virtual ~admVariable()
        {}


    // Member Functions

        // Access

            //- Access time
            inline const admTime& runTime() const;
            inline admTime& runTime();

            //- Access the mesh
            inline const fvMesh& mesh() const;
            
            //- Access variable manager
            inline const admVariableManager& admVars() const;
            inline admVariableManager& admVars();

            //- Variable indices (local and global)
            inline const label& globalIndex() const;
            inline const label& localIndex() const;

            //- Variable type
            inline const varTypeEnum& varType() const;
            
            //- Returns true if this variable is valid (not null)
            inline bool valid() const;

        // Access - semantic shortcuts
        
            //- Variable name
            virtual const word& name() const = 0;

            //- Variable dimensions
            virtual const dimensionSet& dimensions() const = 0;

        // Assign functions

            //- Assign a value at a given cell index
            virtual void assign(const scalar& scalarIn, label cellIndex) = 0;

            //- Assign the dimensions - only checks dimensions
            virtual void assignDims
            (
                const dimensionSet& dimsIn
            ) = 0;
            
            //- Assign a dimensioned value at a given cell index
            virtual void assignDimensioned
            (
                const dimensionedScalar& dScalarIn,
                label cellIndex
            ) = 0;

            //- Assign values to the entire field
            virtual void assignField
            (
                const scalarField& sfIn
            ) = 0;

            //- Assign dimensioned values to the entire field            
            virtual void assignDimensionedField
            (
                const dimensionedScalarField& dsfIn
            ) = 0;

        // Evaluate functions
        //  Evaluate functions for variables return const references where
        //  possible to prevent redundant copying, as opposed to straight up
        //  objects passed in the rest of the framework

            //- Evaluate at a given cell index
            virtual const scalar& evaluate(label cellIndex) const = 0;

            //- Return the name associated with an evaluation
            virtual word evaluateName() const = 0;

            //- Return the dimensions after evaluation
            virtual const dimensionSet& evaluateDims() const = 0;

            //- Evaluate at a given cell index with dimensions
            virtual dimensionedScalar evaluateDimensioned
            (
                label cellIndex
            ) const = 0;

            //- Evaluate the entire internal field
            virtual const scalarField& evaluateField() const = 0;

            //- Evaluate the entire internal field with dimensions
            virtual const dimensionedScalarField&
                evaluateDimensionedField() const = 0;
                
            //- Evaluate the entire Geometric field
            virtual const volScalarField& evaluateGeometricField() const = 0;

            //- True if evaluate is non-zero
            virtual bool evaluateNonZero() const = 0;

        // Ddy functions

            //- Derivative with respect to a given variable at a given cell
            //  index
            virtual scalar ddy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const = 0;

            //- Return the name associated with the derivative with respect to
            //  a given variable
            virtual word ddyName
            (
                const admVariable& wrtVar
            ) const = 0;

            //- Return the dimensions the derivative with respect to a given
            //  variable
            virtual dimensionSet ddyDims
            (
                const admVariable& wrtVar
            ) const = 0;

            //- Derivative with respect to a given variable at a given cell
            //  index with dimensions
            virtual dimensionedScalar ddyDimensioned
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const = 0;

            //- Derivative with respect to a given variable across the entire
            //  internal field
            virtual tmp<scalarField> ddyField
            (
                const admVariable& wrtVar
            ) const = 0;

            //- Derivative with respect to a given variable across the entire
            //  internal field with dimensions
            virtual tmp<dimensionedScalarField> ddyDimensionedField
            (
                const admVariable& wrtVar
            ) const = 0;

            //- True if derivative with respect to a given variable is non-zero
            virtual bool ddyNonZero(const admVariable& wrtVar) const = 0;

        // Housekeeping
            
            //- Correct boundary conditions
            virtual void correctBoundaryConditions() = 0;
            
        // Ancilliary functions

            //- Look up varType number, given a word
            static varTypeEnum wordToVarType(const word&);

            //- Return a word, given a varType
            static word varTypeToWord(const varTypeEnum&);

    // Member operators

        //- Comparison operator - only checks global index
        bool operator==(const admVariable&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admVariableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
