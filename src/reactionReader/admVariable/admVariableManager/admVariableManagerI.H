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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline admTime& admVariableManager::runTime()
{
    return runTime_;
}


inline const admTime& admVariableManager::runTime() const
{
    return runTime_;
}


inline const fvMesh& admVariableManager::mesh() const
{
    return mesh_;
}


inline const IOdictionary& admVariableManager::admVariableDict() const
{
    return admVariableDict_;
}


inline const label& admVariableManager::chainRuleSearchDepth() const
{
    return chainRuleSearchDepth_;
}


inline const scalar& admVariableManager::defaultStandardConvergence() const
{
    return defaultStandardConvergence_;
}


inline const scalar& admVariableManager::defaultImplicitConvergence() const
{
    return defaultImplicitConvergence_;
}


inline const scalar& admVariableManager::defaultStandardMinErrorScale() const
{
    return defaultStandardMinErrorScale_;
}


inline const scalar& admVariableManager::defaultImplicitMinErrorScale() const
{
    return defaultImplicitMinErrorScale_;
}


inline const admVariable& admVariableManager::null() const
{
    return nullStandard_;
}


inline const PtrList<admStandardVariable>& admVariableManager::standard() const
{
    return standard_;
}


inline PtrList<admStandardVariable>& admVariableManager::standard()
{
    return standard_;
}


inline const Foam::admStandardVariable&
    Foam::admVariableManager::standard(label localIndex) const
{
    return standard_[localIndex];
}


inline Foam::admStandardVariable& Foam::admVariableManager::standard
(
    label localIndex
)
{
    return standard_[localIndex];
}


inline const PtrList<admImplicitVariable>& admVariableManager::implicit() const
{
    return implicit_;
}


inline PtrList<admImplicitVariable>& admVariableManager::implicit()
{
    return implicit_;
}


inline const Foam::admImplicitVariable&
    Foam::admVariableManager::implicit(label localIndex) const
{
    return implicit_[localIndex];
}


inline Foam::admImplicitVariable& Foam::admVariableManager::implicit
(
    label localIndex
)
{
    return implicit_[localIndex];
}


inline const PtrList<admDerivedVariable>& admVariableManager::derived() const
{
    return derived_;
}


inline PtrList<admDerivedVariable>& admVariableManager::derived()
{
    return derived_;
}


inline const Foam::admDerivedVariable&
    Foam::admVariableManager::derived(label localIndex) const
{
    return derived_[localIndex];
}


inline Foam::admDerivedVariable& Foam::admVariableManager::derived
(
    label localIndex
)
{
    return derived_[localIndex];
}


inline const UPtrList<admVariable>& admVariableManager::all() const
{
    return all_;
}


inline UPtrList<admVariable>& admVariableManager::all()
{
    return all_;
}


inline const Foam::admVariable& Foam::admVariableManager::all
(
    label globalIndex
) const
{
    return all_[globalIndex];
}


inline Foam::admVariable& Foam::admVariableManager::all(label globalIndex)
{
    return all_[globalIndex];
}


inline UPtrList<admStandardVariable>&
    admVariableManager::changedByUdf()
{
    return changedByUdf_;
}


inline const UPtrList<admStandardVariable>&
    admVariableManager::changedByUdf() const
{
    return changedByUdf_;
}


inline admStandardVariable& admVariableManager::changedByUdf
(
    label index
)
{
    return changedByUdf_[index];
}


inline const admStandardVariable& admVariableManager::changedByUdf
(
    label index
) const
{
    return changedByUdf_[index];
}


inline UPtrList<admStandardVariable>&
    admVariableManager::unchangedByUdf()
{
    return unchangedByUdf_;
}


inline const UPtrList<admStandardVariable>&
    admVariableManager::unchangedByUdf() const
{
    return unchangedByUdf_;
}


inline admStandardVariable& admVariableManager::unchangedByUdf
(
    label index
)
{
    return unchangedByUdf_[index];
}


inline const admStandardVariable& admVariableManager::unchangedByUdf
(
    label index
) const
{
    return unchangedByUdf_[index];
}


inline UPtrList<admImplicitVariable>&
    admVariableManager::implicitAutoSolve()
{
    return implicitAutoSolve_;
}


inline const UPtrList<admImplicitVariable>&
    admVariableManager::implicitAutoSolve() const
{
    return implicitAutoSolve_;
}


inline admImplicitVariable& admVariableManager::implicitAutoSolve
(
    label index
)
{
    return implicitAutoSolve_[index];
}


inline const admImplicitVariable& admVariableManager::implicitAutoSolve
(
    label index
) const
{
    return implicitAutoSolve_[index];
}


inline const IOdictionary& admVariableManager::averagesDict() const
{
    return averagesDict_;
}


inline label admVariableManager::nStandard() const
{
    return standard_.size();
}


inline label admVariableManager::nImplicit() const
{
    return implicit_.size();
}


inline label admVariableManager::nDerived() const
{
    return derived_.size();
}


inline label admVariableManager::nAll() const
{
    return all_.size();
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
