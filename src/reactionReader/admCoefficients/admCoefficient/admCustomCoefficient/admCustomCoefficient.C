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

#include "admCustomCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        admCustomCoefficient,
        0
    );

    addToRunTimeSelectionTable
    (
        admCoefficient,
        admCustomCoefficient,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::admCustomCoefficient::admCustomCoefficient
(
    admTime& runTime,
    const fvMesh& mesh,
    admVariableManager& admVars,
    const dictionary& dict,
    const word& coeffName
)
:
    admCoefficient(admVars, coeffName, dimless),
    runTime_(runTime),
    uniform_
    (
        dict.subDict(coeffName).found("uniform")
      ? Switch(dict.subDict(coeffName).lookup("uniform"))
      : Switch(false)
    ),
    constant_
    (
        dict.subDict(coeffName).found("constant")
      ? Switch(dict.subDict(coeffName).lookup("constant"))
      : Switch(false)
    ),
    constantValueKnown_(false),
    lastField_
    (
        uniform_
      ? NULL
      : new scalarField(mesh.nCells())
    ),
    lastValue_(0),
    lastMilestone_(-1),

    fEqn_
    (
        runTime_.eqns().readEquation
        (
            equation
            (
                dict.subDict(coeffName).lookup("function"),
                coeffName
            )
        )
    ),
    ddtEqn_
    (
        dict.subDict(coeffName).found("ddt")
      ? runTime_.eqns().readEquation
        (
            equation
            (
                dict.subDict(coeffName).lookup("ddt"),
                word("ddt(" + coeffName + ")")
            )
        )
      : label(-1)
    ),
    jEqns_(0),
    jVars_(0),
    dimensionsKnown_(false)
{
    const dictionary& sdict(dict.subDict(coeffName));
    if (sdict.found("Jacobian"))
    {
        dictionary jdict(sdict.subDict("Jacobian"));
        wordList wrtVars(jdict.toc());
        forAll(wrtVars, wrtIndex)
        {
            if (admVars_.found(wrtVars[wrtIndex]))
            {
                label newIndex(jVars_.size());
                jVars_.setSize(newIndex + 1);
                jEqns_.setSize(newIndex + 1);
                jVars_.set
                (
                    newIndex,
                    &admVars_.lookup(wrtVars[wrtIndex])
                );
                
                // read equation and store index
                jEqns_[newIndex] = runTime_.eqns().readEquation
                (
                    equation
                    (
                        jdict.lookup(wrtVars[wrtIndex]),
                        word
                        (
                            "d(" + coeffName + ")|d(" + wrtVars[wrtIndex]
                          + ")"
                        )
                    )
                );
            }
            else
            {
                WarningIn("admCustomCoefficient::admCustomCoefficient")
                    << coeffName << " defines Jacobian with respect to "
                    << "unknown variable " << wrtVars[wrtIndex] << ". "
                    << "Function is being ignored." << endl;
            }
        } // end forAll Jacobian wrtVars
    } // end if Jacobian exists
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dimensionSet& Foam::admCustomCoefficient::dimensions() const
{
    if (!dimensionsKnown_)
    {
        dimensions_.reset
        (
            runTime_.eqns().evaluateDimensions(fEqn_)
        );
        dimensionsKnown_ = true;
    }
    return dimensions_;
}


Foam::scalar Foam::admCustomCoefficient::evaluate(label cellIndex) const
{
    if (constant_)
    {
        if (uniform_)
        {
            if (!constantValueKnown_)
            {
                lastValue_ = runTime_.eqns().evaluateScalar
                (
                    fEqn_, cellIndex
                );
                constantValueKnown_ = true;
            }
            return lastValue_;
        }
        else
        {
            if (!constantValueKnown_)
            {
                return runTime_.eqns().evaluateScalar
                (
                    fEqn_, cellIndex
                );
            }
            return lastField_->operator[](cellIndex);
        }
    }

    // Not constant_
    if (uniform_)
    {
        if (lastMilestone_ != runTime_.milestone())
        {
            lastValue_ = runTime_.eqns().evaluateScalar
            (
                fEqn_, cellIndex
            );
            lastMilestone_ = runTime_.milestone();
        }
        return lastValue_;
    }
    if (lastMilestone_ != runTime_.milestone())
    {
        return runTime_.eqns().evaluateScalar
        (
            fEqn_, cellIndex
        );
    }
    return lastField_->operator[](cellIndex);
}


Foam::tmp<scalarField> Foam::admCustomCoefficient::evaluateField() const
{
    if (constant_)
    {
        if (uniform_)
        {
            if (!constantValueKnown_)
            {
                lastValue_ = runTime_.eqns().evaluateScalar
                (
                    fEqn_, 0
                );
                constantValueKnown_ = true;
            }
            return tmp<scalarField>
            (
                new scalarField
                (
                    mesh_.nCells(),
                    lastValue_
                )
            );
        }
        else
        {
            if (!constantValueKnown_)
            {
                runTime_.eqns().evaluateScalarField
                (
                    * lastField_,
                    fEqn_
                );
                constantValueKnown_ = true;
            }
            return tmp<scalarField>
            (
                new scalarField
                (
                    * lastField_
                )
            );
        }
    }

    // Not constant_
    if (uniform_)
    {
        if (lastMilestone_ != runTime_.milestone())
        {
            lastValue_ = runTime_.eqns().evaluateScalar
            (
                fEqn_, 0
            );
            lastMilestone_ = runTime_.milestone();
        }
        return tmp<scalarField>
        (
            new scalarField
            (
                mesh_.nCells(),
                lastValue_
            )
        );
    }

    if (lastMilestone_ != runTime_.milestone())
    {
        runTime_.eqns().evaluateScalarField
        (
            * lastField_,
            fEqn_
        );
        lastMilestone_ = runTime_.milestone();
    }
    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            * lastField_
        )
    );
    return tReturnMe;
}


Foam::scalar Foam::admCustomCoefficient::ddt(label cellIndex) const
{
    return runTime_.eqns().evaluateScalar(ddtEqn_, cellIndex);
}


Foam::tmp<scalarField> Foam::admCustomCoefficient::ddtField() const
{
    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            mesh_.nCells()
        )
    );
    runTime_.eqns().evaluateScalarField
    (
        tReturnMe(),
        ddtEqn_
    );
    return tReturnMe;
}


bool Foam::admCustomCoefficient::ddtNonZero() const
{
    return (ddtEqn_ >= 0);
}


Foam::scalar Foam::admCustomCoefficient::ddy
(
    const admVariable& wrtVar,
    label cellIndex
) const
{
    forAll(jVars_, jIndex)
    {
        if (jVars_[jIndex].ddyNonZero(wrtVar))
        {
            return scalar
            (
                runTime_.eqns().evaluateScalar
                (
                    jEqns_[jIndex],
                    cellIndex
                ) * jVars_[jIndex].ddy(wrtVar, cellIndex)
            );
        }
    }
    return 0.0;
}


Foam::tmp<scalarField> Foam::admCustomCoefficient::ddyField
(
    const admVariable& wrtVar
) const
{
    tmp<scalarField> tReturnMe
    (
        new scalarField
        (
            mesh_.nCells()
        )
    );
    scalarField& returnMe(tReturnMe());
    
    forAll(jVars_, jIndex)
    {
        if (jVars_[jIndex].ddyNonZero(wrtVar))
        {
            runTime_.eqns().evaluateScalarField
            (
                returnMe,
                jEqns_[jIndex]
            );
            returnMe *= jVars_[jIndex].ddyField(wrtVar);
            return tReturnMe;
        }
    }
    returnMe = 0.0;
    return tReturnMe;
}


bool Foam::admCustomCoefficient::ddyNonZero(const admVariable& wrtVar) const
{
    forAll(jVars_, jIndex)
    {
        if (jVars_[jIndex].ddyNonZero(wrtVar))
        {
            return true;
        }
    }
    return false;
}

// ************************************************************************* //
