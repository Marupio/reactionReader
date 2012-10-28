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
    admReactionReader

Description
    This is a PDE object that OpenFOAM can solve using ODESolver.  It is an
    object-oriented implementation of the Anaerobic Digester Model No. 1 (ADM1)
    as presented in the International Water Association's (IWA) Scientific and
    Technical Report No. 13 (STR13).  The implementation is flexible and can
    accommodate many modifications.

SourceFiles
    admReactionReaderI.H
    admReactionReader.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admReactionReader_H
#define admReactionReader_H

#include "admTime.H"
#include "admVariableManager.H"
#include "admCoefficientManager.H"
#include "admReactionManager.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//template<int matrixSize> class IOadmReactionReaderReferencer;

/*---------------------------------------------------------------------------*\
                            Class admReactionReader Declaration
\*---------------------------------------------------------------------------*/

class admReactionReader
:
    public regIOobject
{
    // Private member data

        //- True when initialize() has been called
        bool initialized_;

protected:

        admTime& runTime_;
        const fvMesh& mesh_;

        //- Dictionary with adm framework-related settings
        IOdictionary admSettingsDict_;
        
        //- Reference to the equation dictionary held by admTime
        const IOdictionary& admEquationDict_;

        //- Managers of model components
        admVariableManager admVars_;
        admCoefficientManager admCoeffs_;
        admReactionManager admReacs_;


        // Output flags
        // These are yes/no settings that change what is output to the console
        // during a run.  They are contained in admSettingsDict.  Other classes
        // can add to this subdictionary.

            //- Reference to outputFlag subdictionary in admSettingsDict,
            //  points to dictionary::null if missing
            const dictionary& outputFlagsDict_;
            
            //- Changes to time and timestep
            bool outputTimeDetails_;

            //- Reaction volumetric averages (default true)
            bool outputReactionAverages_;


        // Non-zero yield reaction lists
        //  These list the reactions that have a non-zero yield for a given
        //  variable.  They reduce the computational effort, especially when
        //  the yield matrix is sparse.

            //- Non-zero yield reactions for all variables:
            //      nonZeroYields_[variable.globalIndex()][i] = reaction
            PtrList<UPtrList<const admReaction> > nonZeroYields_;

            //- Non-zero yield reactions for standard variables:
            //      nonZeroYields_[standardVar.localIndex()][i] = reaction
            PtrList<UPtrList<const admReaction> > nonZeroStandardYields_;

            //- Non-zero yield reactions for implicit variables:
            //      nonZeroYields_[implicitVar.localIndex()][i] = reaction
            PtrList<UPtrList<const admReaction> > nonZeroImplicitYields_;

            //- Non-zero yield reactions for derived variables:
            //      nonZeroYields_[derivedVar.localIndex()][i] = reaction
            PtrList<UPtrList<const admReaction> > nonZeroDerivedYields_;


        // State control
        
            //- Since OpenFOAM doesn't have the ability to return to a previous
            //  timestep, state control is required to enable this behaviour.
            //  The admReactionReader object has no state control variables, but it must
            //  properly manage all the available save state slots.  This is
            //  achieved with a boolList.
            boolList validStates_;

private:

    // Private member functions

        // Initialization - call only once
        
            //- Initialize all component objects, and create internal non-zero
            //  yield lists
            void initialize();
            
        //- Disallow default bitwise copy construct and assignment
        admReactionReader(const admReactionReader&);
        void operator=(const admReactionReader&);
        
        //- Create non-zero index lists for yields
        void createNonZeroYields();

public:

    TypeName("admReactionReader");

    // Constructors

        //- Construct, reading from dictionary
        admReactionReader
        (
            admTime& runTime,
            fvMesh& mesh,
            const word& variableDictName = "admVariableDict",
            const word& coefficientDictName = "admCoefficientDict",
            const word& inhibitionDictName = "admInhibitionDict",
            const word& reactionDictName = "admReactionDict",
            const word& settingsDictName = "admSettingsDict",
            const word& hooksDictName = "admHooksDict"
        );


    // Destructor
    
        virtual ~admReactionReader();


    // Member Functions
    
        // Access

            //- admTime
            inline const admTime& runTime() const;
            inline admTime& runTime();
            
            //- The mesh
            inline const fvMesh& mesh() const;

            //- Dictionary with user-defined settings
            inline const dictionary& admSettingsDict() const;
            inline dictionary& admSettingsDict();

            //- Output flags dictionary (a subdict of admSettingsDict)
            inline const dictionary& outputFlagsDict() const;

            //- Output flag for reaction variable averages
            inline const bool& outputReactionAverages() const;

            //- Dictionary with user-defined equations
            inline const dictionary& admEquationDict() const;

            //- Variable manager            
            inline const admVariableManager& admVars() const;
            inline admVariableManager& admVars();

            //- Coefficient manager
            inline const admCoefficientManager& admCoeffs() const;
            inline admCoefficientManager& admCoeffs();
            
            //- Reaction manager
            inline const admReactionManager& admReacs() const;
            inline admReactionManager& admReacs();
            
            //- Non-zero yields
            inline const PtrList<UPtrList<const admReaction> >&
                nonZeroYields() const;
            inline const PtrList<UPtrList<const admReaction> >&
                nonZeroStandardYields() const;
            inline const PtrList<UPtrList<const admReaction> >&
                nonZeroImplicitYields() const;
            inline const PtrList<UPtrList<const admReaction> >&
                nonZeroDerivedYields() const;


        // Tools
        //  General math functions required by the solver, but useful elsewhere

            //- Calculate the root mean square of the given field
            static scalar calculateRms(const scalarField& sfIn);
            
            //- Calculate the weighted root mean square, given the field of
            //  values and the field of weights
            static scalar calculateRms
            (
                const scalarField& sfIn,
                const scalarField& weights
            );
            
            //- Compare scalarFields a and b, return rms error, scaled by
            //  scalingFactor
            static scalar calculateRmsError
            (
                const scalarField& a,
                const scalarField& b,
                scalar scalingFactor = 1.0
            );

            //- Compare scalarFields a and b using weighted rms error, given
            //  the two fields and the weighting field, scaled by scalingFactor
            static scalar calculateRmsError
            (
                const scalarField& a,
                const scalarField& b,
                const scalarField& weights,
                scalar scalingFactor = 1.0
            );

            //- Returns the "best overall time" based on a comparison between
            //  two given fields, a scaling factor
            static scalar calculateNextDeltaT
            (
                const scalarField& coarseResult,
                const scalarField& fineResult,
                const scalar& convergence,
                const scalarField& weights,
                scalar currentDeltaT,
                scalar scalingFactor = 1.0
            );


        // Housekeeping
        
            //- Sets the time.milestone to force lazy objects to recalculate;
            //  ensures all variable boundary fields are up to date
            void setMilestone();


        // State control
        //  Since OpenFOAM does not easily support undoing a timestep, or
        //  stepping back several timesteps, the save state framework must be
        //  used.
        
            //- Save data to backup variables
            virtual void saveState(const label slot = 0);
            
            //- Erase backup variables
            virtual void clearState(const label slot = 0);
            
            //- Load backup variables into active variables
            virtual void loadState(const label slot = 0);
            
            //- Return number of spots that are occupied
            //  A value of -1 means it doesn't track nStates
            virtual label nStates() const;

            //- Returns true if this is a loadable save state slot
            virtual bool validState(const label slot) const;


        // Required for regIOobject

            //- Write data - do nothing - this IOobject has no input / output
            virtual bool writeData(Ostream&) const;


        // Debugging easy output (output to CSV format)

            //- Report variable values
            void reportVariablesHeader(Ostream& os) const;
            void reportVariables(Ostream& os, const label cellIndex) const;

            //- Report reaction rate values
            void reportReactionRatesHeader(Ostream& os) const;
            void reportReactionRates(Ostream& os, const label cellIndex) const;

            //- Report coefficient values
            void reportCoefficientsHeader(Ostream& os) const;
            void reportCoefficients(Ostream& os, const label cellIndex) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admReactionReaderI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
