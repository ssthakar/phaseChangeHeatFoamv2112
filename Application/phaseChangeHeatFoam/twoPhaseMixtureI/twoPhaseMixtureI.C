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

\*---------------------------------------------------------------------------*/

#include "twoPhaseMixtureI.H"//add (change twoPhaseMixture into twoPhaseMixtureI)
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "dictionary.H"

//#include "IOdictionary.H"
//#include "incompressible/transportModel/transportModel.H"
//#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
//#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void twoPhaseMixtureI::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseMixtureI::twoPhaseMixtureI
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
//    transportModel(U, phi),

    phase1Name_("phase1"),
    phase2Name_("phase2"),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rho1_(nuModel1_->viscosityProperties().lookup("rho")),
    rho2_(nuModel2_->viscosityProperties().lookup("rho")),

    //add ki
    K1_(nuModel1_->viscosityProperties().lookup("K")),
    K2_(nuModel2_->viscosityProperties().lookup("K")),

    //add ci
    C1_(nuModel1_->viscosityProperties().lookup("C")),
    C2_(nuModel2_->viscosityProperties().lookup("C")),

    U_(U),
    phi_(phi),

    alpha1_(U_.db().lookupObject<const volScalarField> (alpha1Name)),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> twoPhaseMixtureI::mu() const
{
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            limitedAlpha1*rho1_*nuModel1_->nu()
          + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
        )
    );
}


tmp<surfaceScalarField> twoPhaseMixtureI::muf() const
{
    surfaceScalarField alpha1f =
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
        )
    );
}


tmp<surfaceScalarField> twoPhaseMixtureI::nuf() const
{
    surfaceScalarField alpha1f =
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
            )/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
        )
    );
}


bool twoPhaseMixtureI::read()
{
    if (transportModel::read())
    {
        if
        (
            nuModel1_().read(subDict(phase1Name_))
         && nuModel2_().read(subDict(phase2Name_))
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
