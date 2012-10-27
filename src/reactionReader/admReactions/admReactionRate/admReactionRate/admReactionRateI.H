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

\*---------------------------------------------------------------------------*/

#include "admReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const admCoefficientManager& admReactionRate::admCoeffs() const
{
    return admCoeffs_;
}


inline admCoefficientManager& admReactionRate::admCoeffs()
{
    return admCoeffs_;
}


inline const admReactionManager& admReactionRate::admReacs() const
{
    return admReacs_;
}


inline admReactionManager& admReactionRate::admReacs()
{
    return admReacs_;
}


inline const UPtrList<const admRateInhibition>
    admReactionRate::rateInhibitions() const
{
    return rateInhibitions_;
}


inline word admReactionRate::inhibitionName() const
{
    return word(name_ + "(rateInhibition)");
}


inline dimensionSet admReactionRate::inhibitionDims() const
{
    return dimless;
}


inline dimensionedScalar admReactionRate::inhibitionDimensioned
(
    label cellIndex
) const
{
    return dimensionedScalar
    (
        inhibitionName(),
        inhibitionDims(),
        inhibition(cellIndex)
    );
}


inline word admReactionRate::inhibitionDdyName
(
    const admVariable& wrtVar
) const
{
    return word("d(" + inhibitionName() + ")|d(" + wrtVar.name() + ")");
}


inline dimensionSet admReactionRate::inhibitionDdyDims
(
    const admVariable& wrtVar
) const
{
    return dimensionSet(inhibitionDims() / wrtVar.dimensions());
}


inline dimensionedScalar admReactionRate::inhibitionDimensioned
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return dimensionedScalar
    (
        inhibitionDdyName(wrtVar),
        inhibitionDdyDims(wrtVar),
        inhibitionDdy(wrtVar, cellIndex)
    );
}


inline word admReactionRate::uninhibitedName() const
{
    return word(name() + "(uninhibited)");
}


inline dimensionSet admReactionRate::uninhibitedDims() const
{
    return dimensions();
}


inline dimensionedScalar admReactionRate::uninhibitedDimensioned
(
    label cellIndex
) const
{
    return dimensionedScalar
    (
        uninhibitedName(),
        uninhibitedDims(),
        uninhibited(cellIndex)
    );
}


inline word admReactionRate::uninhibitedDdyName
(
    const admVariable& wrtVar
) const
{
    return word("d(" + uninhibitedName() + ")|d(" + wrtVar.name() + ")");
}


inline dimensionSet admReactionRate::uninhibitedDdyDims
(
    const admVariable& wrtVar
) const
{
    return uninhibitedDims() / wrtVar.dimensions();
}


inline dimensionedScalar admReactionRate::uninhibitedDdyDimensioned
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    return dimensionedScalar
    (
        uninhibitedDdyName(wrtVar),
        uninhibitedDdyDims(wrtVar),
        uninhibitedDdy(wrtVar, cellIndex)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
