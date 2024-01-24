/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "functionObjectTemplate.H"
#define namespaceFoam  // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NusseltNumberFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    NusseltNumberFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = cccbb51f72797da98bfe2d37f220c5dd4241c6a9
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void NusseltNumber_cccbb51f72797da98bfe2d37f220c5dd4241c6a9(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMesh&
Foam::NusseltNumberFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
NusseltNumberFunctionObject::
NusseltNumberFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
NusseltNumberFunctionObject::
~NusseltNumberFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
NusseltNumberFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read NusseltNumber");
    }

//{{{ begin code
    #line 86 "/home/kazeshini/OpenFOAM/kazeshini-v2112/applications/solvers/phaseChangeHeatFoam_OF1612/Validation/boil/deltaT10/system/controlDict.functions.NusseltNumber"
filePtr.reset(new OFstream("test.dat"));
//}}} end code

    return true;
}


bool
Foam::
NusseltNumberFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute NusseltNumber");
    }

//{{{ begin code
    #line 90 "/home/kazeshini/OpenFOAM/kazeshini-v2112/applications/solvers/phaseChangeHeatFoam_OF1612/Validation/boil/deltaT10/system/controlDict.functions.NusseltNumber"
// declare local constants 
        scalar sigma_ = 0.10; // surface tension between liquid and interface 
        scalar rhog_ = 5; // density of gas
        scalar rhol_ = 200; // density of liquid
        scalar g_ = 9.81;
        scalar Hlg_ = 100e3; // latent heat of vaporization
        scalar kg_ = 1; // thermal conductivity of the gas phase
        scalar mug_ = 0.005;
        scalar DT_ = 5; // overheat
        scalar Cpg_ = 200;

        // critical Taylor wavelenght 
        scalar lambda_ = sqrt( sigma_/((rhol_-rhog_)*g_) ); 
        
				// get curent Temperature Field
        const volScalarField& T
        (
          mesh().lookupObject<volScalarField>("T") 
        );
				
				// compute Berenson's correlation number
				scalar Nub_ = 0.425*pow(((rhog_*(rhol_-rhog_)*g_*Hlg_)/(kg_*mug_*abs(DT_))),0.25)*pow(lambda_,0.75);
				Info << "Berenson's correlation number = " << Nub_ << endl; 
				
				// Compute Klimenko's correlation number
				scalar Beta_ = Cpg_ * DT_ / Hlg_;
				scalar Pr_ = Cpg_ * mug_ / kg_;
				scalar Gr_ = pow(rhog_,2)*g_*pow(lambda_,3)/pow(mug_,2)*(rhol_/rhog_-1);
				scalar Nuk_ = 0.19*pow(Gr_,0.33333333)*pow(Pr_,0.333333)*0.89*pow(Beta_,-0.333333);
				Info << "Klimenko's correlation number = " << Nuk_ << endl; 
				label down = mesh().boundary().findPatchID("down");
				
				// define Numerical Nusselt number
				volScalarField Nusselt
				(
					IOobject
					(
							"Nusselt",
							mesh().time().timeName(),
							mesh(),
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
					),
					mesh(),
					dimensionedScalar("Nusselt", dimless, 0.0)
				);
				Nusselt.boundaryFieldRef()[down] = lambda_/DT_*T.boundaryField()[down].snGrad();
				scalar area = gSum(mesh().magSf().boundaryField()[down]);
				scalar avgNusselt = gSum(Nusselt.boundaryField()[down] * mesh().magSf().boundaryField()[down])/area;
				Info << "Space-averaged Nusselt lambda = " << avgNusselt << "\n" << endl;
				filePtr() << mesh().time().timeName() << " " << avgNusselt << " "<< Nub_ << " "<< Nuk_ << endl;
//}}} end code

    return true;
}


bool
Foam::
NusseltNumberFunctionObject::write()
{
    if (false)
    {
        printMessage("write NusseltNumber");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
NusseltNumberFunctionObject::end()
{
    if (false)
    {
        printMessage("end NusseltNumber");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

