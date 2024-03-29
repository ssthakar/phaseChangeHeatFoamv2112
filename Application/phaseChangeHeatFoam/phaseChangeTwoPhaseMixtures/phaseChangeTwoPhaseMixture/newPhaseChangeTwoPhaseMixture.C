/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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


    description : choose phase change model among existing models

\*---------------------------------------------------------------------------*/

#include "phaseChangeTwoPhaseMixture.H"
#include "twoPhaseMixtureI/twoPhaseMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseChangeTwoPhaseMixture>
Foam::phaseChangeTwoPhaseMixture::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
{
    IOdictionary transportPropertiesDict
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word phaseChangeTwoPhaseMixtureTypeName
    (
        transportPropertiesDict.lookup("phaseChangeTwoPhaseMixture")
    );

    Info<< "Selecting phaseChange model "
        << phaseChangeTwoPhaseMixtureTypeName << endl;


    auto * ctorPtr = componentsConstructorTable(phaseChangeTwoPhaseMixtureTypeName);
    
    
    // componentsConstructorTable::iterator cstrIter =
        // componentsConstructorTablePtr_
            // ->find(phaseChangeTwoPhaseMixtureTypeName);

    if (!ctorPtr)
    {
        FatalErrorIn
        (
            "phaseChangeTwoPhaseMixture::New"
        )   << "Unknown phaseChangeTwoPhaseMixture type "
            << phaseChangeTwoPhaseMixtureTypeName << endl << endl
            << "Valid  phaseChangeTwoPhaseMixtures are : " << endl
            << componentsConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<phaseChangeTwoPhaseMixture>(ctorPtr(U, phi, alpha1Name));
}


// ************************************************************************* //
