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
    Foam::admReaction

Description
    A reaction in the Anaerobic Digester Model (ADM1).  Contains all data
    necessary to describe a single reaction from the admReactionDict file.
    
SourceFiles
    admReactionI.H
    admReaction.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/


#ifndef admReaction_H
#define admReaction_H

#include "admTime.H"
#include "autoPtr.H"
#include "tmp.H"
#include "admVariableManager.H"
#include "admReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class admCoefficientManager;
class admCoefficient;
class admReactionManager;

/*---------------------------------------------------------------------------*\
                           Class admReaction Declaration
\*---------------------------------------------------------------------------*/

class admReaction
{

    // Private data
    
        //- References to other objects
        admTime& runTime_;
        const fvMesh& mesh_;
    
        //- Name of admReaction
        word name_;
        
        //- Reaction rate
        autoPtr<admReactionRate> ratePtr_;
//        admReactionRate& rate_;

        // Yields

            //- Coefficients that give the yield value, indexed by the
            //  variable's globalIndex
            UPtrList<const admCoefficient> yields_;

            //- Sublist - non-zero yields
            labelList nonZeroYields_;

        // Initialization
        
            //- True when initialize() has been performed
            bool initialized_;
            
    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const admReaction&);
        admReaction(const admReaction&);


        //- Read dictionary settings
        void readDict(const dictionary& dict);

public:

    TypeName("admReaction");

    // Constructors

        //- Construct from dictionary
        admReaction
        (
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            admReactionManager& admReacs,
            label standardSize,
            label implicitSize,
            label derivedSize,
            const dictionary& dict,
            const word& reactionName
        );

    // Destructor
    
        virtual ~admReaction()
        {}


    // Member Functions

        // Initialization
        
            //- Initialize - create non-zero yield list
            void initialize();

        // Access

            //- Runtime
            inline const admTime& runTime() const;
            inline admTime& runTime();
            
            //- Mesh
            inline const fvMesh& mesh() const;

            //- Name
            inline const word& name() const;
            
            //- Reaction rate
            inline const admReactionRate& rate() const;
            inline admReactionRate& rate();

            //- Non-zero yield sublist
            inline const labelList& nonZeroYields() const;

        // Yields
        
            //- Return yield coefficient by variable
            inline const admCoefficient& yield(const admVariable& toVar) const;

            //- Return yield coefficient by global variable index
            inline const admCoefficient& yield(label globalIndex) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "admReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
