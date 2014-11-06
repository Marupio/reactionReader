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

Application
    reactionReaderDemo

Description
    A simple application that creates a reactionReader and demonstrates how to
    use it.

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "admReactionReader.H"
#include "HashTableCore.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

    // Note the unusual create time include:
#   include "createAdmTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Create the reaction reader
    admReactionReader model(runTime, mesh);
    
    // Useful references
    admVariableManager& admVars(model.admVars());
    admCoefficientManager& admCoeffs(model.admCoeffs());
    admReactionManager& admReacs(model.admReacs());

    // *** Dealing with variables ***

    // There are three kinds of variables: standard, implicit, derived.
    // You access these through the variable manager.  They are contained in
    // three arrays:
    //
    // * admVars.standard() --> a list of standard variables
    // * admVars.implicit() --> a list of implicit variables
    // * admVars.derived()  --> a list of derived variables
    //
    // There is also an array of all of them together:
    //
    // * admVars.all() --> a list of all variables.
    //
    // The quantities of all of these are readily available:
    Info << "The model has:\n"
        << token::NL << token::TAB
        << admVars.nStandard() << " standard variables;"
        << token::NL << token::TAB
        << admVars.nImplicit() << " implicit variables; and"
        << token::NL << token::TAB
        << admVars.nDerived() << " derived variables,"
        << token::NL
        << "for a total of " << admVars.nAll() << " variables." << endl;
    
    // Table of contents
    // You can get a word list of the variables, but this isn't usually very
    // useful.
    Info << token::NL << "The standard variable names are: "
        << admVars.tocStandard() << endl;
    Info << token::NL << "The implicit variable names are: "
        << admVars.tocImplicit() << endl;
    Info << token::NL << "The derived variable names are: "
        << admVars.tocDerived() << endl;

    // You can search through the variable for one of a specific name:
    if (admVars.found("S_ic"))
    {
        admVariable& tempVar(admVars.lookup("S_ic"));
        Info << "Variable S_ic exists, and has dimensions "
            << tempVar.dimensions() << endl;
    }
    else
    {
        Info << "Variable S_ic doesn't exist." << endl;
    }
    // but it's much faster to use the indexes, see looping, below.
    //
    // There's plenty of search functions, just look through the header for
    // admVariableManager.
    //
    // It's useful to work with variables using a reference.  You can use
    // a general reference (admVariable), but this is limited.  You can't
    // gain direct access to the volScalarField with it.
    admVariable& generalVar(admVars.all(0));

    // Assigning a value (to the internalField):
    generalVar.assign(3.0, 0);
    
    // Retrieving a value:
    Info << "The value in cell 0 is " << generalVar.evaluate(0) << endl;
    
    // There are plenty of other assign, retrieve functions.  See the header
    // for admVariable.
    //
    // Do not create a variable, only create references.  This is a mistake:
    //
    //   admVariable generalVar admVars.all(0);
    //
    // Use this instead:
    //
    //   admVariable& generalVar admVars.all(0);
    //
    // You can do some calculus with the variables:
    admVariable& generalVar2(admVars.all(1));
    Info << "The derivative of " << generalVar.name() << " with respect to "
        << generalVar2.name() << " at cell 0 is = "
        << generalVar.ddy(generalVar2, 0) << endl;
    // This is likely to be zero, unless these are derived variables and they
    // are functions of eachother.

    // It is better to use the specific variable types (e.g.
    // admStandardVariable as opposed to admVariable), you get access to more
    // functions then, including direct access to the volScalarFields:
    admStandardVariable& standardVar(admVars.standard(0));
    volScalarField& vsf(standardVar());
    vsf.internalField()[0] = 4.0;

    // You can loop through variables.  For example:
    forAll(admVars.standard(), varIndex)
    {
        admStandardVariable& sVar(admVars.standard(varIndex));
        Info << "Index " << varIndex << " is " << sVar.name()
            << ", has dimensions " << sVar.dimensions() << endl;
    }

    // *** Dealing with reactions ***

    // The reactionManager has all the reactions.  In fact, it is itself the
    // list, so you can use list functions directly on it for simplicity.
    // For example, the first reaction is given by:
    admReaction& reac0(admReacs[0]);
    
    // And you can loop through the reactions:
    forAll(admReacs, reacIndex)
    {
        admReaction& theReaction(admReacs[reacIndex]);
        Info << "At index " << reacIndex << " it is " << theReaction.name()
            << " and at cell 0 its reaction rate is "
            << theReaction.rate().evaluate(0) << endl;
    }
    
    // With the reaction, you can find out the yields on each variable:
    Info << "At cell 0, the " << reac0.name() << " reaction yields "
        << reac0.yield(generalVar).evaluateDimensioned(0) << " to variable "
        << generalVar.name() << endl;

    // The reaction also gives you the reaction rate:
    admReactionRate& reacRate(reac0.rate());
    // which you can evaluate, take the derivative of, and so on.  See the
    // header file of admReactionRate for more details.

    // *** Dealing with coefficients ***

    // Finally, the coefficients that you defined in admCoefficientsDict are
    // accessed through the coefficientManager.
    Info << "All the coefficients are " << admCoeffs.toc() << endl;
    
    // Coefficients are accessed by name.  They are stored in a hash table, so
    // this is fast.  Use any of the search / find functions from
    // src/OpenFOAM/containers/HashTables/HashTable/HashTable.H, e.g.:
    if (admCoeffs.found("k_a_ac"))
    {
        const admCoefficient& theCoeff(admCoeffs("k_a_ac"));
        Info << "k_a_ac coefficient exists, and at cell 0 its value is "
            << theCoeff.evaluate(0) << endl;
    }
    else
    {
        Info << "k_a_ac coefficient does not exist." << endl;
    }
    // Have a look at the coefficient header file admCoefficient.H.  Notice
    // it inherits admCalculusInterface.  This class defines most of the
    // useful functions you can use, for variables, coefficients, reaction
    // rates, etc..

    // Good luck!
    

    return 0;
}


// ************************************************************************* //
