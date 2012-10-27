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
    Foam::admCalculusInterface

Description
    Abstract base class for fields that perform calculus in the anaerobic
    digester model, ADM1.  This class contains the interface for all objects
    derived from it, including: coefficients, reaction rates, and reaction
    rate inhibitions.
    
SourceFiles
    admCalculusInterface.C
    newadmCalculusInterface.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admCalculusInterface_H
#define admCalculusInterface_H

#include "dimensionedScalar.H"
#include "admTime.H"
#include "fvMesh.H"
#include "admVariableManager.H"
#include "equationVariable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class admCalculusInterface Declaration
\*---------------------------------------------------------------------------*/

class admCalculusInterface
:
    public equationVariable
{

protected:

    //- Reference to runTime
    admTime& runTime_;

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Reference to admVariables for Jacobian
    admVariableManager& admVars_;

    //- Name of calculus field
    word name_;

    //- Dimensions of calculus field
    mutable dimensionSet dimensions_;

public:

    //- Runtime type information
    TypeName("admCalculusInterface");


    // Constructors

        //- Construct from components
        admCalculusInterface
        (
            admVariableManager& admVars,
            const word& name,
            const dimensionSet& dimensions
        );

    // Destructor

        virtual ~admCalculusInterface()
        {}


    // Member functions

        // Access

            //- runTime
            inline admTime& runTime()
            {
                return runTime_;
            }
            
            //- runTime const
            inline const admTime& runTime() const
            {
                return runTime_;
            }

            //- mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to admVariables
            admVariableManager& admVars()
            {
                return admVars_;
            }

            //- Reference to admVariables (const)
            const admVariableManager& admVars() const
            {
                return admVars_;
            }

            //- Return the name of the calculus field
            virtual const word& name() const
            {
                return name_;
            }

            //- Return the dimensions of the calculus field
            virtual const dimensionSet& dimensions() const
            {
                return dimensions_;
            }

        // Evaluation functions

            //- Evaluate a given cell index
            virtual scalar evaluate(label cellIndex) const = 0;

            //- Return the name associated with an evaluation
            virtual const word& evaluateName() const
            {
                return name();
            }

            //- Return the dimensions after evaluation
            virtual const dimensionSet& evaluateDims() const
            {
                return dimensions();
            }

            //- Evaluate at a given cell index with dimensions
            inline dimensionedScalar evaluateDimensioned
            (
                label cellIndex
            ) const
            {
                return dimensionedScalar
                (
                    evaluateName(),
                    evaluateDims(),
                    evaluate(cellIndex)
                );
            }

            //- Evaluate the entire internal field
            virtual tmp<scalarField> evaluateField() const = 0;

            //- Evaluate the entire internal field with dimensions
            tmp<dimensionedScalarField> evaluateDimensionedField() const;

            //- True if evaluate is non-zero
            virtual bool evaluateNonZero() const
            {
                return true;
            }

        // ddy functions

            //- Derivative with respect to a given variable at a given cell
            //  index
            virtual scalar ddy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const
            {
                return 0.0;
            }

            //- Return the name associated with the derivative with respect to
            //  a given variable
            virtual word ddyName
            (
                const admVariable& wrtVar
            ) const
            {
                return word("d(" + name() + ")|d(" + wrtVar.name() + ")");
            }

            //- Return the dimensions the derivative with respect to a given
            //  variable
            virtual dimensionSet ddyDims
            (
                const admVariable& wrtVar
            ) const
            {
                return dimensions() / wrtVar.dimensions();
            }

            //- Derivative with respect to a given variable at a given cell
            //  index with dimensions
            inline dimensionedScalar ddyDimensioned
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const
            {
                return dimensionedScalar
                (
                    ddyName(wrtVar),
                    ddyDims(wrtVar),
                    ddy(wrtVar, cellIndex)
                );
            }

            //- Derivative with respect to a given variable across the entire
            //  internal field
            virtual tmp<scalarField> ddyField
            (
                const admVariable& wrtVar
            ) const;

            //- Derivative with respect to a given variable across the entire
            //  internal field with dimensions
            tmp<dimensionedScalarField> ddyDimensionedField
            (
                const admVariable& wrtVar
            ) const;

            //- True if derivative with respect to a given variable is non-zero
            virtual bool ddyNonZero(const admVariable& wrtVar) const
            {
                return false;
            }

        // Required by equationVariable class
            
            virtual label lookupComponentIndex(const word) const
            {
                return 0;
            }
            
            virtual scalar evaluateScalar
            (
                const label componentIndex,
                label cellIndex,
                const label geoIndex
            ) const
            {
                return evaluate(cellIndex);
            }

            virtual void evaluateScalarField
            (
                scalarField& result,
                const label componentIndex,
                const label geoIndex
            ) const
            {
                result = evaluateField();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
