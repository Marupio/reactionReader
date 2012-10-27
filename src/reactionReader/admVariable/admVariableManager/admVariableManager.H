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
    Foam::admVariableManager

Description
    Holds all the anaerobic digester variables, including: standard, implicit,
    and derived.  Provides interface functions for these.

SourceFiles
    admVariableManagerI.H
    admVariableManager.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admVariableManager_H
#define admVariableManager_H

#include "admStandardVariable.H"
#include "admImplicitVariable.H"
#include "admDerivedVariable.H"
#include "IOdictionary.H"
#include "fvCFD.H"
#include "objectRegistry.H"
//#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class admReactionReader;

/*---------------------------------------------------------------------------*\
                       Class admVariableManager Declaration
\*---------------------------------------------------------------------------*/

class admVariableManager
{

    // Private member data

        //- Quick reference
        admReactionReader& model_;
        admTime& runTime_;
        const fvMesh& mesh_;

        // User-defined settings
        
            //- Dictionary with user-defined settings
            IOdictionary admVariableDict_;

            //- Chain-rule search depth
            label chainRuleSearchDepth_;

        // Default settings
        
            //- Convergence
            scalar defaultStandardConvergence_;
            scalar defaultImplicitConvergence_;
            
            //- Minimum error scaling value
            scalar defaultStandardMinErrorScale_;
            scalar defaultImplicitMinErrorScale_;

        // Variable lists

            //- Variable definitions for returning null
            admStandardVariable nullStandard_;
            admImplicitVariable nullImplicit_;
            admDerivedVariable nullDerived_;

            //- Variables
            PtrList<admStandardVariable> standard_;
            PtrList<admImplicitVariable> implicit_;
            PtrList<admDerivedVariable> derived_;
            
            //- Full variable index
            UPtrList<admVariable> all_;


        // Variable list subsets

            //- Index list of standard variables marked with changedByUdf
            UPtrList<admStandardVariable> changedByUdf_;

            //- Index list of standard variables not marked with changedByUdf
            UPtrList<admStandardVariable> unchangedByUdf_;
            
            //- Index list of implicit variables marked with autoSolve
            UPtrList<admImplicitVariable> implicitAutoSolve_;


        // Input / output
        
            //- Dictionary holding the volume-weighted averages for each
            //  variable. Only updated before runTime.write()
            mutable IOdictionary averagesDict_;


        // Initialization
        
            //- true when initialize() has already been called
            bool initialized_;
        

        // State control

            //- Standard variable restore point fields (saves all .old fields)
            PtrList<PtrList<volScalarField> > standardRestore_;

            //- Implicit variable restore point fields (saves all .old fields)
            PtrList<PtrList<volScalarField> > implicitRestore_;
            
            //- State control of derived variables is handled by the individual
            //  objects.


    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        admVariableManager(const admVariableManager&);
        void operator=(const admVariableManager&);

        //- Read admVariableDict
        void readDict();

public:

    TypeName("admVariableManager");


    // Constructors

        //- Construct, reading from dictionary
        admVariableManager
        (
            admReactionReader& model,
            const word& variableDictName = "admVariableDict"
        );


    // Destructor
    
        virtual ~admVariableManager()
        {}


    // Member Functions

        // Initialization - call only once
            //- Perform any initializiation routines
            void initialize();


        // Access

            // Quick reference
            
                //- Runtime
                inline admTime& runTime();
                inline const admTime& runTime() const;

                //- Mesh
                inline const fvMesh& mesh() const;


            // User-defined settings
            
                //- Dictionary with user-defined settings
                inline const IOdictionary& admVariableDict() const;

                //- Chain-rule search depth
                inline const label& chainRuleSearchDepth() const;

            // Default settings

                //- Convergence criteria
                inline const scalar& defaultStandardConvergence() const;
                inline const scalar& defaultImplicitConvergence() const;
                
                //- Minimum error scale factor
                inline const scalar& defaultStandardMinErrorScale() const;
                inline const scalar& defaultImplicitMinErrorScale() const;


            // Variable lists

                //- Return null variable
                inline const admVariable& null() const;

                //- Standard variables
                inline const PtrList<admStandardVariable>& standard() const;
                inline PtrList<admStandardVariable>& standard();

                //- Standard variables by local index
                inline const admStandardVariable&
                    standard(label localIndex) const;
                inline admStandardVariable& standard(label localIndex);
                
                //- Implicit variables
                inline const PtrList<admImplicitVariable>& implicit() const;
                inline PtrList<admImplicitVariable>& implicit();

                //- Implicit variable by local index
                inline const admImplicitVariable&
                    implicit(label localIndex) const;
                inline admImplicitVariable& implicit(label localIndex);

                //- Derived variables
                inline const PtrList<admDerivedVariable>& derived() const;
                inline PtrList<admDerivedVariable>& derived();
                
                //- Derived variable by local index
                inline const admDerivedVariable&
                    derived(label localIndex) const;
                inline admDerivedVariable& derived(label localIndex);
                
                //- All variables
                inline const UPtrList<admVariable>& all() const;
                inline UPtrList<admVariable>& all();
                
                //- Access any variable by globalIndex
                inline const admVariable& all(label globalIndex) const;
                inline admVariable& all(label globalIndex);


            // Variable list subsets

                //- Index list of standard variables marked with changedByUdf
                inline UPtrList<admStandardVariable>&
                    changedByUdf();
                inline const UPtrList<admStandardVariable>&
                    changedByUdf() const;
                
                //- changedByUdf sublist by index
                inline admStandardVariable& changedByUdf
                (
                    label index
                );
                inline const admStandardVariable& changedByUdf
                (
                    label index
                ) const;

                //- Index list of standard variables not marked with
                //  changedByUdf
                inline UPtrList<admStandardVariable>&
                    unchangedByUdf();
                inline const UPtrList<admStandardVariable>&
                    unchangedByUdf() const;
                
                //- unchangedByUdf sublist by index
                inline admStandardVariable& unchangedByUdf
                (
                    label index
                );
                inline const admStandardVariable& unchangedByUdf
                (
                    label index
                ) const;

                //- Index list of implicit variables marked with autoSolve
                inline UPtrList<admImplicitVariable>&
                    implicitAutoSolve();
                inline const UPtrList<admImplicitVariable>&
                    implicitAutoSolve() const;
                
                //- implicitAutoSolve by index
                inline admImplicitVariable& implicitAutoSolve
                (
                    label index
                );
                inline const admImplicitVariable& implicitAutoSolve
                (
                    label index
                ) const;

            //- Dictionary holding the volume-weighted averages for each
            //  variable. Only updated before runTime.write()
            inline const IOdictionary& averagesDict() const;


            // Number of variables
                
                //- Return number of standard variables
                inline label nStandard() const;
                
                //- Return number of implicit variables
                inline label nImplicit() const;
                
                //- Return number of derived variables
                inline label nDerived() const;
                
                //- Return number of all variables
                inline label nAll() const;


        // Searching
        
            //- True if varName is a standard variable
            bool foundStandard(const word& varName) const;

            //- True if varName is an implicit variable
            bool foundImplicit(const word& varName) const;

            //- True if varName is a derived variable
            bool foundDerived(const word& varName) const;

            //- True if varName is any variable
            bool found(const word& varName) const;

            //- Find and return the standard variable
            const admStandardVariable& lookupStandard
            (
                const word& varName
            ) const;
            admStandardVariable& lookupStandard(const word& varName);

            //- Find and return the implicit variable
            const admImplicitVariable& lookupImplicit
            (
                const word& varName
            ) const;
            admImplicitVariable& lookupImplicit(const word& varName);

            //- Find and return the derived variable
            const admDerivedVariable& lookupDerived(const word& varName) const;
            admDerivedVariable& lookupDerived(const word& varName);

            //- Find and return the variable from any list
            const admVariable& lookup(const word& varName) const;
            admVariable& lookup(const word& varName);
            
            //- List available variables
            wordList tocStandard() const;
            wordList tocImplicit() const;
            wordList tocDerived() const;
            wordList toc() const;


        // Variable limits

            //- Apply limits to a geometric internal field, while maintaining
            //  a weighted average throughout the internal field by removing
            //  or adding an offset value to the rest of the values
            void massConservingApplyLimits
            (
                scalarField& sf,
                const scalarField& weights,
                const scalar& upperLimit,
                const scalar& lowerLimit,
                const word& varName = word("unknown")
            );
        
            //- Apply limits to a standard variable across the entire mesh,
            //  ensuring mass is neither created nor destroyed
            void massConservingApplyLimits(admStandardVariable& var);

            //- Apply limits to an implicit variable across the entire mesh,
            //  ensuring mass is neither created nor destroyed
            void massConservingApplyLimits(admImplicitVariable& var);

            //- Apply limits to all standard variables across the entire mesh,
            //  ensuring mass is neither created nor destroyed
            void massConservingApplyStandardLimits();

            //- Apply limits to all implicit variables across the entire mesh,
            //  ensuring mass is neither created nor destroyed
            void massConservingApplyImplicitLimits();

            //- Apply limits to any type of variable across the internal field
            void applyLimits(admVariable& var);

            //- Apply limits to a standard variable across the internal field
            void applyLimits(admStandardVariable& var);

            //- Apply limits to an implicit variable across the internal field
            void applyLimits(admImplicitVariable& var);

            //- Apply limits to any type of variable variable at a cell index
            void applyLimits
            (
                admVariable& var,
                label cellIndex
            );

            //- Apply limits to an standard variable at a cell index
            void applyLimits
            (
                admStandardVariable& var,
                label cellIndex
            );

            //- Apply limits to an implicit variable at a cell index
            void applyLimits
            (
                admImplicitVariable& var,
                label cellIndex
            );

            //- Apply limits to the standard variables across internal field
            void applyStandardLimits();

            //- Apply limits to the standard variables at specified cell index
            void applyStandardLimits(label cellIndex);

            //- Apply limits to the implicit variables across internal field
            void applyImplicitLimits();

            //- Apply limits to the implicit variables at specified cell index
            void applyImplicitLimits(label cellIndex);

            //- Apply limits to the derived variables across internal field
            void applyDerivedLimits();

            //- Apply limits to the derived variables at specified cell index
            void applyDerivedLimits(label cellIndex);


        // Internal field savespots
        //  Apparently some algorithms need the internal fields of variables
        //  to be saved for temporary purposes
        
            //- Save standard variable internal fields
            void saveStandardInternalFields
            (
                PtrList<scalarField>& saveSpot
            );

            //- Save implicit variable internal fields
            void saveImplicitInternalFields
            (
                PtrList<scalarField>& saveSpot
            );

            //- Save derived variable internal fields
            void saveDerivedInternalFields
            (
                PtrList<scalarField>& saveSpot
            );

            //- Save all variable internal fields
            void saveAllInternalFields
            (
                PtrList<scalarField>& saveSpot
            );


        // Housekeeping

            //- Correct the averages of each variable to the time directory
            void updateAverages() const;
            
            //- Ensure all the boundary conditions are up-to-date
            void correctStandardBoundaryConditions();
            void correctImplicitBoundaryConditions();
            void correctDerivedBoundaryConditions();
            void correctAllBoundaryConditions();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admVariableManagerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
