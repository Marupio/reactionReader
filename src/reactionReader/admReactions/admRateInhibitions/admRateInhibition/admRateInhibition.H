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
    Foam::admRateInhibition

Description
    Abstract base class for inhibition that acts as a multiplier to the
    reaction rates (rate) in the anaerobic digester model, ADM1.  This class
    contains the interface for all inhibition multiplier rates.
    
SourceFiles
    admRateInhibition.C
    newAdmRateInhibition.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admRateInhibition_H
#define admRateInhibition_H

#include "refCount.H"
#include "autoPtr.H"
#include "dimensionedScalar.H"
#include "admVariableManager.H"
#include "admCoefficientManager.H"
#include "admReactionJumpTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class admRateInhibition Declaration
\*---------------------------------------------------------------------------*/

class admRateInhibition
:
    public admCalculusInterface
{

protected:

    //- Reference to coefficients manager
    admCoefficientManager& admCoeffs_;
    
    //- The full field values go here
    mutable scalarField lastValue_;

    // Milestone - if different from admTime's milestone, evaluation
    //  requires a new calculation
    mutable unsigned int lastMilestone_;

private:
    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const admRateInhibition&);
        admRateInhibition(const admRateInhibition&);

public:

    //- Runtime type information
    TypeName("admRateInhibition");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            admRateInhibition,
            dictionary,
            (
                admTime& runTime,
                admVariableManager& admVars,
                admCoefficientManager& admCoeffs,
                const dictionary& dict,
                const word& inhibitionName
            ),
            (runTime, admVars, admCoeffs, dict, inhibitionName)
        );

    // Constructors

        //- Construct given components
        admRateInhibition
        (
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            const word& name
        );

    // Selectors

        //- Return a reference to a new admRateInhibition
        //  *** Located in newAdmRateInhibition.C ***
        static autoPtr<admRateInhibition> New
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            const dictionary& dict,
            const word& inhibitionName
        );

    // Destructor

        virtual ~admRateInhibition()
        {}

    // Member functions

        // Access
        
            //- Return the coefficients manager
            inline const admCoefficientManager& admCoeffs() const;
            inline admCoefficientManager& admCoeffs();

        // Evaluation functions
        //  These functions automatically provide lazy-evaluation.  It is
        //  recommended you do not define these in your derived classes,
        //  rather, define the inhibition functions, on which the evaluation
        //  functions rely for their values.

            //- Evaluate a given cell index
            virtual scalar evaluate(label cellIndex) const;

            //- Evaluate the entire internal field
            virtual tmp<scalarField> evaluateField() const;

            //- True if evaluate is non-zero
            virtual bool evaluateNonZero() const;

        // Inhibition functions
        //  These are regular evaluate functions, without lazy-evaluation.
        //  Define these in derived classes.

            //- Evaluate inhibition at a given cell index
            virtual scalar inhibition(label cellIndex) const = 0;

            //- Evaluate inhibition throughout entire internal field
            virtual tmp<scalarField> inhibitionField() const = 0;

            //- True if inhibition is non-zero
            virtual bool inhibitionNonZero() const;

        // Ddy functions

            virtual scalar ddy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const = 0;

            virtual tmp<scalarField> ddyField
            (
                const admVariable& wrtVar
            ) const = 0;

            virtual bool ddyNonZero(const admVariable& wrtVar) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admRateInhibitionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
