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
    Foam::admTemperatureDependentExponentialRatioCoefficient

Description
    A coefficient for the Anaerobic Digester Model (ADM1) that varies with
    temperature according to the relation:

    \f[
        C = A * exp
            \left (
                B * \left ( \frac{1}{T_{base}} - \frac{1}{T} \right )
            \right )
    \f]

    where:
        A and B are numerical constants;
        \f$ T_{base} \f$ is a base temperature; and
        T is the current temperature.
    
SourceFiles
    admTemperatureDependentExponentialRatioCoefficient.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admTemperatureDependentExponentialRatioCoefficient_H
#define admTemperatureDependentExponentialRatioCoefficient_H

#include "admCoefficient.H"
// #include "volFieldsFwd.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class admTemperatureDependentExponentialRatioCoefficient Declaration
\*---------------------------------------------------------------------------*/

class admTemperatureDependentExponentialRatioCoefficient
:
    public admCoefficient
{
    // Private data

        //- The full field values go here
        mutable scalarField lastValue_;

        // Milestone - if different from admTime's milestone, evaluation
        //  requires a new calculation
        mutable unsigned int lastMilestone_;
        
        //- Numerical constants
        scalar A_;
        scalar B_;
        scalar T_base_;

        //- Associated variables
        const admVariable& varT_;
        
        //- If varT_ is derived, this is not null, affects ddt calculations
        admDerivedVariable const * derivedVarTPtr_;


public:

    //- Runtime type information
    TypeName("temperatureDependentExponentialRatio");


    // Constructors

        //- Construct from dictionary, given name
        admTemperatureDependentExponentialRatioCoefficient
        (
            admTime& runTime,
            const fvMesh& mesh,
            admVariableManager& admVars,
            const dictionary& dict,
            const word& coeffName
        );

    // Destructor

        virtual ~admTemperatureDependentExponentialRatioCoefficient()
        {}


    // Member functions
    
        // Access

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
