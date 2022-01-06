/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixture::phaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
//------------------22.01.06----------------------//
    const word& alpha1Name
//------------------------------------------------//
)
:
    incompressibleTwoPhaseMixture(U, phi),
    phaseChangeTwoPhaseMixtureCoeffs_(optionalSubDict(type + "Coeffs")),
    pSat_("pSat", dimPressure, lookup("pSat")),
//------------------------------------add 21.12.27-----------------------------------//
    TSat_("TSat", dimensionSet (0, 0, 0, 1, 0, 0, 0), lookup("TSat")),
    TSatLocal_(readBool(lookup("TSatLocal"))),
    Hfg_("Hfg", dimensionSet (0, 2, -2, 0, 0, 0, 0), lookup("Hfg")),
    R_("R", dimensionSet (0, 2, -2, -1, 0, 0, 0), lookup("R"))

//-----------------------------------------------------------------------------------//

{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1()*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}

//------------------------------------add 21.12.27-----------------------------------//

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotT() const
{
	   volScalarField rhoCp =  rho1()*C1()*alpha1_ + rho2()*C2()*(1.0-alpha1_); //2phase rhoCp 정의
	   //volScalarField TCoeff = (Hfg()+(C1()-C2())*TSat_)/rhoCp;
	   volScalarField TCoeff = Hfg()/rhoCp;										//2phase에 의한 source term 계수
	   //volScalarField TCoeff = Hfg()*(alpha1_/C1() + (1.0 - alpha1_)/C2());
	   Pair<tmp<volScalarField> > mDotT = this->mDotT();

	    return Pair<tmp<volScalarField>> (TCoeff*mDotT[0], TCoeff*mDotT[1]);
}

Foam::volScalarField
Foam::phaseChangeTwoPhaseMixture::TSatLocal() const
{
	if (TSatLocal_)
	{
		//Info <<"TSatlocal" << endl;
	    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
	    volScalarField oneByT = 1.0/TSat() - R_/Hfg_*log(max(p/pSat_,1E-08));

	    return volScalarField
	    (
	           1.0/oneByT
	    );
	}
	else
	{
		//Info <<"TSat" << endl;
		volScalarField one
		(
			IOobject
			(
				"one",
				alpha1_.mesh(),
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			alpha1_.mesh(),
			dimensionedScalar ("one",dimless, 1.0)
	    );

		return volScalarField
		(
			one*TSat_
	    );
	}
}

//-----------------------------------------------------------------------------------//

void Foam::phaseChangeTwoPhaseMixture::correct()
{
    incompressibleTwoPhaseMixture::correct();
}


bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");
        lookup("pSat") >> pSat_;
//------------------------------------add 21.12.27-----------------------------------//
        lookup("TSat") >> TSat_;
        lookup("TSatLocal") >> TSatLocal_;
        lookup("Hfg") >> Hfg_;
        lookup("R") >> R_;

//-----------------------------------------------------------------------------------//
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
