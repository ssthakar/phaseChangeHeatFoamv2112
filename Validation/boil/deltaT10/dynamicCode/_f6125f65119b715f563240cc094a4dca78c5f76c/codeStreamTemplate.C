/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 23 "/home/kazeshini/OpenFOAM/kazeshini-v2112/applications/solvers/phaseChangeHeatFoam_OF1612/Validation/boil/deltaT10/0/T.#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_f6125f65119b715f563240cc094a4dca78c5f76c(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 39 "/home/kazeshini/OpenFOAM/kazeshini-v2112/applications/solvers/phaseChangeHeatFoam_OF1612/Validation/boil/deltaT10/0/T.#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>(dict);
		const fvMesh& mesh = refCast<const fvMesh>(d.db());
		scalar lambda = 0.0786844;
    scalar deltaT = 5.;
    scalar Twall = 505;
    scalar Tliquid = 500.;
    scalarField T(mesh.nCells(), Tliquid);
		forAll(T, i)
		{
			const scalar x = mesh.C()[i][0];
			const scalar y = mesh.C()[i][1];
			if ( y <= lambda/128*(4+cos(2*constant::mathematical::pi*x/lambda)) )
			{
				T[i] = -deltaT*y / (lambda/128 * (4 + cos(2*constant::mathematical::pi*x/lambda))) + Twall;
			}
		}
		T.writeEntry("", os);
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

