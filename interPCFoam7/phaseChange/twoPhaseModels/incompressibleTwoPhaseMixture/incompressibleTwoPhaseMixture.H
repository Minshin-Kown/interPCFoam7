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
    Foam::incompressibleTwoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    incompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseMixture_H
#define incompressibleTwoPhaseMixture_H

#include "../../transportModels/transportModel/transportModel.H"
#include "../../transportModels/viscosityModels/viscosityModel/viscosityModel.H"
#include "../../twoPhaseModels/twoPhaseMixture/twoPhaseMixture/twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class incompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseMixture
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

//---------------------add for upgrade 21.12.17-------------------//
    	  word phase1Name_;
	      word phase2Name_;
//----------------------------------------------------------------//        
        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

//---------------------add for upgrade 21.12.17-------------------//

        dimensionedScalar K1_;		//thermal conductivity
	dimensionedScalar K2_;

        dimensionedScalar C1_;		//specific heat
	dimensionedScalar C2_;

//----------------------------------------------------------------//

        const volVectorField& U_;
        const surfaceScalarField& phi_;

//---------------------add for upgrade 21.12.17-------------------//
	const volScalarField& alpha1_;	//twoPhaseMixture에 정의
//----------------------------------------------------------------//

        volScalarField nu_;


    // Protected Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    TypeName("incompressibleTwoPhaseMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
//---------------------add for upgrade 21.12.17-------------------//
	    ,const word& alpha1Name = "alpha1" 	//twoPhaseMixture에 정의
//----------------------------------------------------------------//
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseMixture()
    {}


    // Member Functions

//---------------------add for upgrade 22.01.03-------------------//
        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }
//----------------------------------------------------------------//

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        }

//---------------------add for upgrade 21.12.17-------------------//

        
        //- Return const-access to phase1 specific heat
        const dimensionedScalar& C1() const //add ci
        {
            return C1_;
        }

        //- Return const-access to phase2 specific heat
        const dimensionedScalar& C2() const
        {
            return C2_;
        }
        //- Return const-access to phase1 heat conductivity
        const dimensionedScalar& K1() const //add ki
        {
        	return K1_;
        }

        //- Return const-access to phase2 heat conductivity
        const dimensionedScalar& K2() const
        {
        	return K2_;
        }
      
//----------------------------------------------------------------//

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
