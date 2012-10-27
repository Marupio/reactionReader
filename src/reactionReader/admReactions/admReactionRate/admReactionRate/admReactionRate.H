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
    Foam::admReactionRate

Description
    Abstract base class for reaction rates in the anaerobic digester model,
    ADM1. This class contains the interface for all reaction rates.
    
SourceFiles
    admReactionRate.C
    newAdmReactionRate.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admReactionRate_H
#define admReactionRate_H

#include "autoPtr.H"
#include "tmp.H"
#include "dimensionedScalar.H"
#include "admTime.H"
#include "fvMesh.H"
#include "admVariableManager.H"
#include "admCoefficientManager.H"
#include "admCalculusInterface.H"
#include "admRateInhibition.H"
#include "admReactionJumpTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class admReactionManager;

/*---------------------------------------------------------------------------*\
                        Class admReactionRate Declaration
\*---------------------------------------------------------------------------*/

class admReactionRate
:
    public admCalculusInterface
{

protected:

    // Protected member data
    
        //- Reference to coefficient manager for constants in dictionaries
        admCoefficientManager& admCoeffs_;

        //- Reference to reaction manager required for inhibitions
        admReactionManager& admReacs_;
        
        //- Associated rateInhibition(s), if any
        UPtrList<const admRateInhibition> rateInhibitions_;
        
        // Non-zero index lists

            //- inhibitionDdy - calculated using chainrule.  Non-zero indices
            //  correspond to an entire product group that is non-zero:
            //      I0 I1 ... dIi / dy ... In-1 In
            //  where the given index is i (the rateInhibition in the
            //  derivative):
            //      nonZeroInhibitionDdy[wrtVar][i] = rateInhibitionIndex
            labelListList nonZeroInhibitionDdy_;
            
            //- True when non-zero index lists have been initialized
            bool initialized_;
        
        //- The full field values go here
        mutable scalarField lastEvaluateValue_;
        mutable scalarField lastInhibitionValue_;

        // Milestone - if different from admTime's milestone, evaluation
        //  requires a new calculation
        mutable unsigned int lastEvaluateMilestone_;
        mutable unsigned int lastInhibitionMilestone_;

private:

    // Private member functions
    
        //- Read dictionary entries (inhibition list)
        void readDict(const dictionary& dict, const word& reactionName);

        //- Disallow default bitwise copy construct and assignment
        void operator=(const admReactionRate&);
        admReactionRate(const admReactionRate&);
        
public:

    //- Runtime type information
    TypeName("admReactionRate");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            admReactionRate,
            dictionary,
            (
                admTime& runTime,
                admVariableManager& admVars,
                admCoefficientManager& admCoeffs,
                admReactionManager& admReacs,
                const dictionary& dict,
                const word& reactionName
            ),
            (runTime, admVars, admCoeffs, admReacs, dict, reactionName)
        );

    // Constructors

        //- Construct given components
        admReactionRate
        (
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            const dictionary& dict,
            const word& name
        );

        //- Construct copy
//        admReactionRate(const admReactionRate&);

    // Selectors

        //- Return a reference to a new admReactionRate
        //  *** Located in newAdmReactionRate.C ***
        static autoPtr<admReactionRate> New
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            const dictionary& dict,
            const word& reactionName
        );

    // Destructor

        virtual ~admReactionRate()
        {}

    // Member functions

        // Access

            //- Return the coefficient manager
            inline const admCoefficientManager& admCoeffs() const;

            inline admCoefficientManager& admCoeffs();
            
            //- Return the reaction manager
            inline const admReactionManager& admReacs() const;

            inline admReactionManager& admReacs();

            //- Return associated rate inhibitions
            inline const UPtrList<const admRateInhibition>
                rateInhibitions() const;

        // Initialization
        
            //- Initialize the object (fills non-zero index lists)
            void initialize();

        // Inhibition functions
        //  All the values of reactionRateInhibition multiplied together

            //- Inhibition at a given cell index
            virtual scalar inhibition(label cellIndex) const;

            //- Return the name of the inhibition
            inline word inhibitionName() const;

            //- Return the dimensions of the inhibition (dimless)
            inline dimensionSet inhibitionDims() const;

            //- Inhibition at a given cell index with dimensions
            inline dimensionedScalar inhibitionDimensioned
            (
                label cellIndex
            ) const;
            
            //- Inhibition across the entire field
            virtual tmp<scalarField> inhibitionField() const;

            //- Inhibition across the entire field with dimensions
            tmp<dimensionedScalarField> inhibitionDimensionedField() const;

            //- True if inhibition is non-zero
            virtual bool inhibitionNonZero() const;

        // Inhibition ddy functions
        //  The derivative of the inhibition with respect to a given variable

            //- The inhibitionDdy at a given cell index
            virtual scalar inhibitionDdy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            //- Return the name associated with the inhibitionDdy
            inline word inhibitionDdyName(const admVariable& wrtVar) const;

            //- Return the dimensions of the inhibitionDdy
            inline dimensionSet inhibitionDdyDims
            (
                const admVariable& wrtVar
            ) const;

            //- The inhibitionDdy at a given cell index with dimensions
            inline dimensionedScalar inhibitionDimensioned
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            
            //- The inhibitionDdy across the entire field
            virtual tmp<scalarField> inhibitionDdyField
            (
                const admVariable& wrtVar
            ) const;

            //- The inhibitionDdy across the entire field with dimensions
            tmp<dimensionedScalarField> inhibitionDdyDimensionedField
            (
                const admVariable& wrtVar
            ) const;

            //- True if the inhibitionDdy is non-zero
            virtual bool inhibitionDdyNonZero(const admVariable& wrtVar) const;

        // Uninhibited value functions
        //  The reactionRate, not including the effects of inhibition

            //- Uninhibited value at a given cell index
            virtual scalar uninhibited(label cellIndex) const = 0;

            //- Return the name of the uninhibited value
            inline word uninhibitedName() const;

            //- Return the dimensions of the uninhibited value
            inline dimensionSet uninhibitedDims() const;

            //- Uninhibited value at a given cell index with dimensions
            inline dimensionedScalar uninhibitedDimensioned
            (
                label cellIndex
            ) const;
            
            //- Uninhibited value across the entire field
            virtual tmp<scalarField> uninhibitedField() const = 0;

            //- Uninhibited value across the entire field with dimensions
            tmp<dimensionedScalarField>
                uninhibitedDimensionedField() const;

            //- True if uninhibited value is non-zero
            virtual bool uninhibitedNonZero() const
            {
                return true;
            }

        // Uninhibited ddy functions
        //  The reactionRate, not including the effects of inhibition

            //- Uninhibited ddy at a given cell index
            virtual scalar uninhibitedDdy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const = 0;

            //- Return the name of the uninhibited ddy
            inline word uninhibitedDdyName(const admVariable& wrtVar) const;

            //- Return the dimensions of the uninhibited ddy
            inline dimensionSet uninhibitedDdyDims
            (
                const admVariable& wrtVar
            ) const;

            //- Uninhibited ddy at a given cell index with dimensions
            inline dimensionedScalar uninhibitedDdyDimensioned
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;
            
            //- Uninhibited ddy across the entire field
            virtual tmp<scalarField> uninhibitedDdyField
            (
                const admVariable& wrtVar
            ) const = 0;

            //- Uninhibited ddy across the entire field with dimensions
            tmp<dimensionedScalarField> uninhibitedDdyDimensionedField
            (
                const admVariable& wrtVar
            ) const;

            //- True if uninhibited ddy is non-zero
            virtual bool uninhibitedDdyNonZero
            (
                const admVariable& wrtVar
            ) const = 0;


        // Evaluation functions
        //  The reactionRate, including the effects of inhibition

            //- Evaluate a given cell index
            virtual scalar evaluate(label cellIndex) const;

            //- Evaluate the entire internal field
            virtual tmp<scalarField> evaluateField() const;

            //- True if evaluate is non-zero
            virtual bool evaluateNonZero() const;

        // ddy functions
        //  The derivative of the reactionRate with respect to a given
        //  variable, including the effects of inhibition

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

#include "admReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
