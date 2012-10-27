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
    Foam::admReactionManager

Description
    Manages all the reactions for the Anaerobic Digester Model (ADM) defined in
    the constant/admReactionDict file.

SourceFiles
    admReactionManagerI.H
    admReactionManager.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admReactionManager_H
#define admReactionManager_H

#include "admTime.H"
#include "admRateInhibition.H"
#include "admCoefficientManager.H"
#include "admReaction.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class admReactionManager Declaration
\*---------------------------------------------------------------------------*/

class admReactionManager
:
    public PtrList<admReaction>
{

// Forward declarations
//class admReaction;

    // Private member data

        //- References to other objects
        admTime& runTime_;
        admVariableManager& admVars_;
        admCoefficientManager& admCoeffs_;
        
        //- Dictionary with the rate inhibition definitions
        IOdictionary admInhibitionDict_;

        //- Dictionary with the reaction definitions
        IOdictionary admReactionDict_;

        //- Multiplier rate inhibitions
        HashPtrTable<const admRateInhibition> inhibitions_;

        //- True when initialize() has been performed
        bool initialized_;
        
    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        admReactionManager(const admReactionManager&);
        void operator=(const admReactionManager&);
        
        //- Read dictionary
        void readDict();

public:

    TypeName("admReactionManager");

    // Constructors

        //- Construct, reading from dictionary
        admReactionManager
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            const word& inhibitionDictName = "admInhibitionDict",
            const word& reactionDictName = "admReactionDict"
        );

    // Destructor
    
        virtual ~admReactionManager()
        {}

    // Member Functions

        // Access

            //- Dictionary with rate inhibition definitions
            inline const IOdictionary& admInhibitionDict() const
            {
                return admInhibitionDict_;
            }

            //- Dictionary with reaction definitions
            inline const IOdictionary& admReactionDict() const
            {
                return admReactionDict_;
            }

            //- Rate inhibition pointer list
            const HashPtrTable<const admRateInhibition>& inhibitions() const
            {
                return inhibitions_;
            }

        // Initialization
        void initialize();

        // Rate inhibition functions

            //- Return a the full list of rate inhibition names
            wordList inhibitionToc() const;

            //- True if inhibitionName exists
            bool inhibitionFound
            (
                const word& inhibitionName
            ) const;

            //- Return the inhibition by name            
            const admRateInhibition& inhibitionLookup
            (
                const word& inhibitionName
            ) const;

        // Reactions
        
            //- Return a list of all reaction names
            wordList toc() const;

            //- True if reactionName exists
            bool found(const word& reactionName) const;

            //- Access a reaction by name
            const admReaction& lookup(const word& reactionName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "admReactionManagerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
