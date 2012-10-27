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

#include "admReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const admTime& admReaction::runTime() const
{
    return runTime_;
}


inline admTime& admReaction::runTime()
{
    return runTime_;
}


inline const fvMesh& admReaction::mesh() const
{
    return mesh_;
}


inline const word& admReaction::name() const
{
    return name_;
}


inline const admReactionRate& admReaction::rate() const
{
    return ratePtr_();
}


inline admReactionRate& admReaction::rate()
{
    return ratePtr_();
}


inline const labelList& admReaction::nonZeroYields() const
{
    return nonZeroYields_;
}


inline const admCoefficient& admReaction::yield(const admVariable& toVar) const
{
    return yields_[toVar.globalIndex()];
}


inline const admCoefficient& admReaction::yield(label globalIndex) const
{
    return yields_[globalIndex];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
