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

#include "admDerivedVariable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const bool& admDerivedVariable::dimensionsKnown() const
{
    return dimensionsKnown_;
}


inline const Switch& admDerivedVariable::suppressOutput() const
{
    return suppressOutput_;
}


inline const label& admDerivedVariable::funEqn() const
{
    return funEqn_;
}


inline const label& admDerivedVariable::ddtEqn() const
{
    return ddtEqn_;
}


inline const UPtrList<admVariable>& admDerivedVariable::rVars() const
{
    return rVars_;
}


inline const wordList& admDerivedVariable::rNames() const
{
    return rNames_;
}


inline const labelList& admDerivedVariable::rEqns() const
{
    return rEqns_;
}


inline const UPtrList<admVariable>& admDerivedVariable::jVars() const
{
    return jVars_;
}


inline const wordList& admDerivedVariable::jNames() const
{
    return jNames_;
}


inline const labelList& admDerivedVariable::jEqns() const
{
    return jEqns_;
}


inline const UPtrList<admVariable>& admDerivedVariable::cVars() const
{
    return cVars_;
}


inline const labelListList& admDerivedVariable::cEqns() const
{
    return cEqns_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
