/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::phaseChangeTwoPhaseMixture

Description

SourceFiles
    phaseChangeTwoPhaseMixture.C
    phaseChangeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeTwoPhaseMixture_H
#define phaseChangeTwoPhaseMixture_H

//#include "immiscibleIncompressibleTwoPhaseMixture.H"	//for devide interface and phangeChangeTwoPhaseMixture
#include "../../twoPhaseModels/incompressibleTwoPhaseMixture/incompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class phaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeTwoPhaseMixture
:
    public incompressibleTwoPhaseMixture	//changed from immiscibleIncompressibleTwoPhaseMixture
{
protected:

    // Protected data

        dictionary phaseChangeTwoPhaseMixtureCoeffs_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;

//----------------------add 21.12.27--------------------//

        //- Saturation vapour Temperature
        dimensionedScalar TSat_; //add
        //- true if Saturation vapour Temperature based on pressure
        bool TSatLocal_; //add
        //- latent Heat
        dimensionedScalar Hfg_;//add
        //- gas specific constant
        dimensionedScalar R_;//add

//-------------------------------------------------------//



public:

    //- Runtime type information
    TypeName("phaseChangeTwoPhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeTwoPhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
//------------------------22.01.05---------------------------//		
		const word& alpha1Name
//-------------------------------5---------------------------//
            ),
            (U, phi, alpha1Name)    //----22.01.05
        );


    // Constructors

        //- Construct from components
        phaseChangeTwoPhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
//------------------------22.01.05---------------------------//
	    const word& alpha1Name = "alpha1"
//-----------------------------------------------------------//
        );

        //- Disallow default bitwise copy construction
        phaseChangeTwoPhaseMixture(const phaseChangeTwoPhaseMixture&);


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeTwoPhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
//------------------------22.01.05---------------------------//
            const word& alpha1Name = "alpha1"
//----------------------------------------------------------//
        );


    //- Destructor
    virtual ~phaseChangeTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

//------------------------------------------add 21.12.27-------------------------------------//

	//- Return const-access to the saturation vapour Temperature
        const dimensionedScalar& TSat() const //add
        {
        	return	TSat_;
        }

        //- Return const-access to the latent heat
        const dimensionedScalar& Hfg() const //add
        {
        	return Hfg_;
        }

//-------------------------------------------------------------------------------------------//
        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const = 0;

//------------------------------------------add 21.12.27-------------------------------------//

	virtual Pair<tmp<volScalarField> > mDotT() const = 0;
        
        Pair<tmp<volScalarField> > vDotT() const;
        
        //- Return const-access to the saturation vapour Temperature
        volScalarField TSatLocal() const;

//-------------------------------------------------------------------------------------------//

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply (p - pSat)
        Pair<tmp<volScalarField>> vDotP() const;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeTwoPhaseMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
