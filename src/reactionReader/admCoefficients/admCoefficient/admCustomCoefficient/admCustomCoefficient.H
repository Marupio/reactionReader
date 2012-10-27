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
    Foam::admCustomCoefficient

Description
    A coefficient for the Anaerobic Digester Model (ADM1) with a custom
    spatial / temporal relation.  The user supplies the function, its time
    derivative (ddt), and its derivative w.r.t. all the admVariables
    (Jacobian).  Furthermore, the user specifies if the coefficient varies with
    space.
    
SourceFiles
    admCustomCoefficient.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admCustomCoefficient_H
#define admCustomCoefficient_H

#include "admCoefficient.H"
// #include "volFieldsFwd.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class admCustomCoefficient Declaration
\*---------------------------------------------------------------------------*/

class admCustomCoefficient
:
    public admCoefficient
{
    // Private data

        //- Reference to adm runTime, gives access to equationReader
        admTime& runTime_;

        //- When true, this coefficient does not vary spatially
        const Switch uniform_;
        
        //- When true, this coefficient does not vary with time
        const Switch constant_;
        
        //- True when constant, and the value has been determined
        mutable bool constantValueKnown_;

        //- The latest full field values go here
        mutable scalarField * lastField_;
        
        //- The latest single value, if uniform is true, goes here
        mutable scalar lastValue_;

        // Milestone - if different from admTime's milestone, evaluation
        //  requires a new calculation
        mutable unsigned int lastMilestone_;

        //- Main equation index
        label fEqn_;
        
        //- Ddt equation index
        label ddtEqn_;

        //- Jacobian equation indices
        labelList jEqns_;

        //- Jacobian wrt variables
        UPtrList<const admVariable> jVars_;
        
        //- Dimensions are only known after initialization
        mutable bool dimensionsKnown_;

public:

    //- Runtime type information
    TypeName("custom");


    // Constructors

        //- Construct from dictionary, given name
        admCustomCoefficient
        (
            admTime& runTime,
            const fvMesh& mesh,
            admVariableManager& admVars,
            const dictionary& dict,
            const word& coeffName
        );

    // Destructor

        virtual ~admCustomCoefficient()
        {}

    // Member functions

        // Access

        // Access

            //- Returns true if the coefficient does not change spatially
            virtual bool uniform() const
            {
                return uniform_;
            }
            
            //- Returns true if the coefficient is invariant with all variables
            virtual bool constant() const
            {
                return constant_;
            }

            //- Return the dimensions of the coefficient
            virtual const dimensionSet& dimensions() const;

        // Evaluation

            //- Evaluate a given cell index
            virtual scalar evaluate(label cellIndex) const;

            //- Evaluate the entire internal field
            virtual tmp<scalarField> evaluateField() const;

        // ddt evaluation

            //- Evaluate ddt at a given cell index
            virtual scalar ddt(label cellIndex) const;

            //- Evaluate ddt over the entire internal field
            virtual tmp<scalarField> ddtField() const;

            //- True if ddt is non-zero
            virtual bool ddtNonZero() const;

        // ddy functions

            //- Derivative with respect to a given variable at a given cell
            //  index
            virtual scalar ddy
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

            //- True if derivative with respect to a given variable is non-zero
            virtual bool ddyNonZero(const admVariable& wrtVar) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
