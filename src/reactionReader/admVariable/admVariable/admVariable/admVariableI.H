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

#include "admVariable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const admTime& admVariable::runTime() const
{
    return runTime_;
}


inline admTime& admVariable::runTime()
{
    return runTime_;
}


inline const fvMesh& admVariable::mesh() const
{
    return mesh_;
}


inline const admVariableManager& admVariable::admVars() const
{
    return admVars_;
}


inline admVariableManager& admVariable::admVars()
{
    return admVars_;
}


inline const label& admVariable::globalIndex() const
{
    return globalIndex_;
}


inline const label& admVariable::localIndex() const
{
    return localIndex_;
}


inline const admVariable::varTypeEnum& admVariable::varType() const
{
    return varType_;
}


inline bool admVariable::valid() const
{
    return (globalIndex_ >= 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
