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
    Foam::admTime

Description
   Manages an anaerobic digester simulation.  Includes an equationReader for
   yield functions and custom models, such as admCustomReactionRate.

SourceFiles
    admTime.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admTime_H
#define admTime_H

#include "Time.H"
#include "IOEquationReader.H"
#include "fvCFD.H"
#include "dimensionedScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class admTime Declaration
\*---------------------------------------------------------------------------*/

class admTime
:
    public Time
{
    // Private data

        //- equationReader for universal equation handling
        IOEquationReader eqns_;

        //- Dictionary with any additional equations
        IOdictionary admEquationDict_;
        
        //- True after doubleDeltaT was called, false after halveDeltaT.
        //  Initializes to false.
        bool doubled_;
        
        //- Timestep limits
        scalar maxDeltaT_;
        scalar minDeltaT_;

        //- Holds the last time value when output was generated, used for
        //  minimum spaced output
        mutable scalar lastOutputTime_;
        
        //- Minimum spacing between output (in run time seconds)
        scalar minimumOutputSpacing_;

        //- Milestone counter
        //  The milestone counter is used to enable the ADM framework's lazy
        //  evaluation.  When a calculation is requested, an ADM object will
        //  check the milestone to see if its number differs from the last time
        //  the calculation was performed.
        unsigned int milestone_;

        //- State control
        scalarList timeValue_a;
        labelList timeIndex_a;
        scalarList deltaT_a;
        scalarList deltaTSave_a;
        scalarList deltaT0_a;
        boolList deltaTchanged_a;
        labelList outputTimeIndex_a;
        boolList outputTime_a;
        boolList doubled_a;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        admTime(const admTime&);

        //- Disallow default bitwise assignment
        void operator=(const admTime&);

public:

    //- Runtime type information
    TypeName("admTime");

    // Constructors

        //- Construct from objectRegistry arguments
        admTime
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const fileName& systemName = "system",
            const fileName& constantName = "constant",
            const word& equationsDictName = "admEquationDict"
        );

    // Destructor

        virtual ~admTime()
        {}


    // Member Functions

        // Access

            //- Const access to the equationReader
            const IOEquationReader& eqns() const
            {
                return eqns_;
            }

            //- Access to the equationReader
            IOEquationReader& eqns()
            {
                return eqns_;
            }
            
            const IOdictionary& eqnDict() const
            {
                return admEquationDict_;
            }

            const bool& doubled() const
            {
                return doubled_;
            }

        // Extensions to equationReader
            
            //- Read a scalar or an equation.  Accepts any of these formats:
            //      keyword     scalar;
            //      keyword     "equation";
            //      keyword     [dimless] "equation";
            //      keyword     ignoredWord [dimless] "equation";
            //  Evaluates immediately, and discards the equation.
            scalar readScalarOrEquation
            (
                Istream& is,
                const word& varName
            );
            
            //- Read a scalar, dimensionedScalar, nameless dimensioned scalar
            //  Accepts any of these formats:
            //      keyword     scalar;
            //      keyword     [dimensionSet] scalar;
            //      keyword     ignoredWord [dimensionSet] scalar;
            //      keyword     "equation";
            //      keyword     [dimensionSet] "equation";
            //      keyword     ignoredWord [dimensionSet] "equation";
            //  Evaluates immediately and discards the equation
            dimensionedScalar readDimensionedScalarOrEquation
            (
                Istream& is,
                const word& varName
            );

        // Milestone functions
        
            //- Set milestone - increment the counter by one
            void setMilestone()
            {
                milestone_++;
            }
            
            //- Return latest milestone value
            const unsigned int& milestone() const
            {
                return milestone_;
            }

        // Edit

            //- Double the time step - only used to check the accuracy of the
            //  half-value time step.  Does not check if it exceeds maxDeltaT.
            //  Only for use with adaptive timestepping.
            void doubleDeltaT();

            //- Return to a regular time step
            void halveDeltaT();

            //- Set the timestep before the coarse step
            void setDeltaTCoarse(const scalar deltaT);

            //- Reset time step
            virtual void setDeltaT(const scalar);
            
            //- Reset the time, time-index, deltaT and deltaT0, no checks on
            //  deltaT.  Use to return to a previously valid state.
            virtual void setTime
            (
                const scalar& newTimeValue,
                const label& newTimeIndex,
                const scalar& newDeltaT,
                const scalar& newDeltaTSave,
                const scalar& newDeltaT0,
                const bool& newDeltaTchanged,
                const label& newOutputTimeIndex,
                const bool& newOutputTime
            );

        // Database functions

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType
            ) const;

        // State control
        
            //- Save data to backup variables
            void saveState(const label slot = 0);
            
            //- Erase backup variables
            void clearState(const label slot = 0);
            
            //- Load backup variables into active variables
            void loadState(const label slot = 0);
            
            //- Return number of spots that are occupied
            //  A value of -1 means it doesn't track nStates
            label nStates() const;

            //- Returns true if this is a loadable save state slot
            bool validState(const label slot) const;

    // Member operators

        //- Fake operator++ - does the same thing a operator++ except it does
        //  not allow output.  Used for timestep-doubling adaptive timestep.
        admTime& plusPlusNoOutput();

        //- Prefix increment
        virtual admTime& operator++();

        //- Postfix increment
        virtual admTime& operator++(int);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
