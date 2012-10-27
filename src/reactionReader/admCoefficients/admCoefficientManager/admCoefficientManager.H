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
    Foam::admCoefficientManager

Description
    Holds all the anaerobic digester coefficients.  These are the model
    constants, and includes almost  everything in tables 2.2, 2.3, and 2.4 in
    the IWA Scientific and Technical Report Number 13, such as:
    
    @verbatim
        C_i         Carbon content of component i;
        N_i         Nitrogen content of component i;
        f_pro_su    Catabolic yield of product on substrate
        H_gas       Gas law constant
        K_a_acid    Acid-base equilibrium coefficients
        K_H         Henry's law coefficient
        k_dec       Decay coefficient
        Y_substrate Yield of biomass on substrate
    @endverbatim
    
    It does not include reaction rates (rate_j), inhibition multipliers (I), or
    yields in the Peterson Matrix (nu_ij).  Uses an abstract base class
    "admCoefficient", for which derived entities exist, including constant,
    temperatureDependentExponentialRatio (allows for temperature variation),
    and custom (allows for anything else, such as inhibitions that act on the
    coefficients).

SourceFiles
    admCoefficientManagerI.H
    admCoefficientManager.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admCoefficientManager_H
#define admCoefficientManager_H

#include "admTime.H"
#include "admCoefficient.H"
#include "admZeroCoefficient.H"
#include "admOneCoefficient.H"
#include "IOdictionary.H"
#include "fvCFD.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class admCoefficientManager Declaration
\*---------------------------------------------------------------------------*/

class admCoefficientManager
:
    public HashPtrTable<admCoefficient>
{

    // Private member data

        admTime& runTime_;
        const fvMesh& mesh_;
        admVariableManager& admVars_;
        
        //- Dictionary with user-defined settings
        IOdictionary admCoefficientDict_;
        
        //- For functions returning zero or one (these are dimless)
        admZeroCoefficient zero_;
        admOneCoefficient one_;
        
    // Private member functions
    
        //- Disallow default bitwise copy construct and assignment
        admCoefficientManager(const admCoefficientManager&);
        void operator=(const admCoefficientManager&);

        //- Read admCoefficientDict
        void readDict();

public:

    TypeName("admCoefficientManager");

    // Constructors

        //- Construct, reading from dictionary
        admCoefficientManager
        (
            admTime& runTime,
            const fvMesh& mesh,
            admVariableManager& admVars,
            const word& coefficientDictName = "admCoefficientDict"
        );

    // Destructor
    
        virtual ~admCoefficientManager()
        {}

    // Member Functions

        // Access

            //- Dictionary with user-defined settings
            inline const IOdictionary& admCoefficientDict() const
            {
                return admCoefficientDict_;
            }
            
            //- Return the zero coefficient (const)
            inline const admCoefficient& zero() const
            {
                return zero_;
            }

            //- Return the zero coefficient
            inline admCoefficient& zero()
            {
                return zero_;
            }

            //- Return the one coefficient (const)
            inline const admCoefficient& one() const
            {
                return one_;
            }

            //- Return the one coefficient
            inline admCoefficient& one()
            {
                return one_;
            }

        // Searching
            
            //- Returns a list of the coefficient names
            wordList toc() const;

            //- Find coefficient by name
            const admCoefficient& operator()(const word& coeffName) const;

        // Add new coefficients
        //  These functions take a dictionary as input.  Acceptable formats
        //  include:
        //      lookupName  word; // readCoefficient only
        //      lookupName  scalarOrEquation;
        //      lookupName  [dimensionSet] scalarOrEquation;
        //      lookupName  ignoredWord [dimensionSet] scalarOrEquation;
        //      lookupName  {sub-dictionary}
        //  The first format returns the named coefficient.  The remaining
        //  formats create a new coefficient with the name newCoeffName,
        //  and return that.  If no name is provided for newCoeffName,
        //  lookupName is the new coefficient name.

            //- Add a new coefficient using a dictionary format.  If coeffName
            //  is given, the newly created coefficient is given this name.
            const admCoefficient& addNew
            (
                const dictionary& dict,
                const word& lookupName,
                const word& newCoeffName = word::null
            );
            
            //- Read a coefficient from dictionary
            const admCoefficient& readCoefficient
            (
                const dictionary& dict,
                const word& lookupName,
                const word& newCoeffName = word::null
            );

/*
        // State control
        
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
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "admCoefficientManagerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
