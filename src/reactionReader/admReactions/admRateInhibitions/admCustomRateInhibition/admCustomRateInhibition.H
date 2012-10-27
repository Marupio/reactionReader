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
    Foam::admCustomRateInhibition

Description
    A reaction rate inhibition multiplier for the Anaerobic Digester Model
    (ADM1) given by a user-defined function and its Jacobian.
    
SourceFiles
    admCustomRateInhibition.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef admCustomRateInhibition_H
#define admCustomRateInhibition_H

#include "admRateInhibition.H"
//#include "admCustomCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class admCustomRateInhibition Declaration
\*---------------------------------------------------------------------------*/

class admCustomRateInhibition
:
    public admRateInhibition
{
    // Private data

        //- A custom coefficient has all the functionality we need - just have
        //  to interface with it
        const admCoefficient& inhibitionValue_;

public:

    //- Runtime type information
    TypeName("custom");


    // Constructors

        //- Construct from dictionary, given name
        admCustomRateInhibition
        (
            admTime& runTime,
            admVariableManager& admVars,
            admCoefficientManager& admCoeffs,
            const dictionary& dict,
            const word& inhibitionName
        );

    // Destructor

        virtual ~admCustomRateInhibition()
        {}


    // Member functions
    
        // Inhibition functions

            //- Evaluate inhibition at a given cell index
            virtual scalar inhibition(label cellIndex) const;

            //- Evaluate inhibition throughout entire internal field
            virtual tmp<scalarField> inhibitionField() const;

            //- True if inhibition is non-zero
            virtual bool inhibitionNonZero() const;

        // Ddy functions

            virtual scalar ddy
            (
                const admVariable& wrtVar,
                label cellIndex
            ) const;

            virtual tmp<scalarField> ddyField
            (
                const admVariable& wrtVar
            ) const;

            virtual bool ddyNonZero(const admVariable& wrtVar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
