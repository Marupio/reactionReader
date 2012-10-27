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
    Foam::admDerivedVariable

Description
    Anaerobic digester derived variable - a variable that can be expressed as
    a single function of standard and implicit variables.

SourceFiles
    admDerivedVariableI.H
    admDerivedVariable.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admDerivedVariable_H
#define admDerivedVariable_H

#include "admVariable.H"
#include "equation.H"
#include "equationVariable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class admDerivedVariable Declaration
\*---------------------------------------------------------------------------*/

class admDerivedVariable
:
    public admVariable,
    public equationVariable
{

    // Private member data

        //- When true, IOobject is set to NO_WRITE
        Switch suppressOutput_;

        //- Pointer to the volScalarField containing the data
        mutable autoPtr<volScalarField> dataPtr_;

        //- Dimensions cannot be learned until all variables are read in and
        //  the function is evaluated for the first time
        bool dimensionsKnown_;

        //- Function equation index
        label funEqn_;
        
        //- Function time derivative equation index; -1 = no ddt exists
        label ddtEqn_;

        // Reverse-function-related
        
            //- Reverse function dictionary pointer - temporary during object
            //  construction
            dictionary * rdict_;

            //- Reverse function associated variables
            UPtrList<admVariable> rVars_;
            
            //- Reverse function names
            wordList rNames_;
            
            //- Reverse function equation indices
            labelList rEqns_;

        // Jacobian-related (ddy)

            //- Jacobian dictionary pointer - temporary during object
            //  construction
            dictionary * jdict_;

            //- Jacobian "with respect to" variable list
            UPtrList<admVariable> jVars_;
            
            //- Jacobian equation names
            wordList jNames_;
        
            //- Jacobian equation indices
            labelList jEqns_;
        
        // Chain-rule-related
        //  Derived variables may miss a non-zero ddy() evaluation unless they
        //  know to what variables they are related through multiple
        //  applications of the chain rule
        
            //- cVars_[i] - the ith chained variable
            UPtrList<admVariable> cVars_;

            //- cEqns_[i][j] - gives a series of equations (indexed by j) that
            //  must be evaluated to get the total Jacobian of the ith variable
            labelListList cEqns_;

        // Milestone - if different from admTime's milestone, evaluation
        //  requires a new calculation
        mutable unsigned int lastMilestone_;
        

private:

    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        admDerivedVariable(const admDerivedVariable&);
        void operator=(const admDerivedVariable&);

    // Indices

        //- Read jdict's of variables - removed from constructor because
        //  admVars must be fully populated before jdict can be read
        void createJacobianIndices();

        //- Read rdict's of variables - returns true if reverse functions
        //  are defined for this variable
        void createReverseFunctionIndices();

        //- Return list of chain-rule-linked variables (recursive)
        //  searchDepth is input
        //  varTypes, varIndices, are indices for variables that are linked via
        //      the chain rule.  These are output values.
        //  eqnIndices is the associated list of Jacobian equations in the
        //      equationReader.  This is also an output value.
        //  operations is the number of operations required to calculate the
        //      Jacobian - used to defeat circular dependency
        void findChainedVariables
        (
            const label& searchDepth,
            UPtrList<admVariable>& varList,
            labelListList& eqnIndices,
            labelList& operations
        ) const;

        // Used by findChainedVariables - appends the 'next' lists to the end
        // of the 'current' lists.  Duplicate entries (same varType, same
        // varIndex) are overwritten by the one with the fewest operations.
        void appendChainedLists
        (
            UPtrList<admVariable> nextVars,
            labelListList& nextEqnIndices,
            labelList& nextOperations,
            UPtrList<admVariable> currentVars,
            labelListList& currentEqnIndices,
            labelList& currentOperations
        ) const;

public:

    TypeName("admDerivedVariable");

    // Constructors

        //- Construct null - internal use only - to permit returning null
        admDerivedVariable
        (
            admTime& runTime,
            admVariableManager& admVars
        );

        //- Construct from dictionary
        admDerivedVariable
        (
            admTime& runTime,
            admVariableManager& admVars,
            const dictionary& dict,
            const label globalIndex,
            const label localIndex,
            const word& name
        );

    // Destructor

        virtual ~admDerivedVariable()
        {}


    // Member Functions

        // Access

            //- Return the dimensionsKnown flag
            inline const bool& dimensionsKnown() const;

            //- Return the suppressOutput flag
            inline const Switch& suppressOutput() const;

            //- Return function equation index
            inline const label& funEqn() const;
        
            //- Return time derivative equation index; -1 = no ddt exists
            inline const label& ddtEqn() const;

            //- Return reverse function associated variables
            inline const UPtrList<admVariable>& rVars() const;
            
            //- Return reverse function names
            inline const wordList& rNames() const;
            
            //- Return reverse function equation indices
            inline const labelList& rEqns() const;

            //- Return Jacobian "with respect to" variable list
            inline const UPtrList<admVariable>& jVars() const;

            //- Return Jacobian equation names
            inline const wordList& jNames() const;
        
            //- Return Jacobian equation indices
            inline const labelList& jEqns() const;
        
            //- Return chain variable list:
            //      cVars()[i] - the ith chained variable
            inline const UPtrList<admVariable>& cVars() const;

            //- Return the chain rule equation indices:
            //      cEqns()[i][j] - gives a series of equations (indexed by j)
            //      that must be evaluated to get the total Jacobian of the
            //      ith variable
            inline const labelListList& cEqns() const;

        // Access - semantic shortcuts
        
            //- Variable name
            virtual const word& name() const;

            //- Variable dimensions
            virtual const dimensionSet& dimensions() const;

        // Housekeeping

            //- Correct boundary conditions
            virtual void correctBoundaryConditions();
        
            //- Create indices for Jacobian and reverse functions - shouldn't
            //  be called twice, but safe if it is
            void createIndices();
        
            //- Initialize - cannot be initialized until all variables are
            //  known, variable cannot be used until initialized - shouldn't
            //  be called twice, but safe if it is
            void initialize();
            
            //- Create chain rule indices - this starts the search for this
            //  variable's chain rule variables
            void createChainRuleIndices(const label chainRuleSearchDepth);

        // Assign functions
        //  admDerivedVariable behaviour: assign the value, and run the reverse
        //  functions

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

        // Ddt functions
        //  Only available for derived variables.  Non-derived variables have
        //  non-zero ddt's, but there is no clear way to determine this.  It is
        //  up to the solver

            //- Evaluate ddt at a given cell index
            scalar ddt(label cellIndex) const;

            //- Return the name associated with a ddt evaluation
            word ddtName() const;

            //- Return the dimensions after a ddt evaluation
            dimensionSet ddtDims() const;

            //- Evaluate ddt at a given cell index with dimensions
            dimensionedScalar ddtDimensioned
            (
                label cellIndex
            ) const;

            //- Evaluate ddt over the entire internal field
            tmp<scalarField> ddtField() const;

            //- Evaluate ddt over the entire internal field with dimensions
            tmp<dimensionedScalarField> ddtDimensionedField() const;

            //- True if ddt is non-zero
            bool ddtNonZero() const;

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

            //- Access the data directly - const only for derivedVariable
            //  This does not update the variable.  It only gives access to the
            //  last calculated value
            inline const volScalarField& operator()() const;

        // Required by equationVariable interface

            virtual label lookupComponentIndex(const word) const;
            
            virtual scalar evaluateScalar
            (
                const label componentIndex,
                label cellIndex,
                const label geoIndex
            ) const;

            virtual void evaluateScalarField
            (
                scalarField& result,
                const label componentIndex,
                const label geoIndex
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admDerivedVariableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
