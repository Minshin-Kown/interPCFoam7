Info<< "--------------------------------\n" << endl;
Info<< "Create fluid fields\n" << endl;


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field alpha1\n" << endl;
volScalarField alpha1
(
    IOobject
    (
        "alpha1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Adding to Kinetic energy field KE\n" << endl;
/*
volScalarField KE
(
 	"KE",
	0.5*magSqr(U)
);
*/

volScalarField KE
(
     IOobject
     (
      	"KE",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
     ),
     0.5*magSqr(U)
);     


#include "createPhi.H"

//-------------------------------------------------------------------//

Info<< "Reading fluid thermophysical properties\n" << endl;


autoPtr<fluidThermo> myfluidThermoPtr (fluidThermo::New(mesh));

fluidThermo& FluidThermo = myfluidThermoPtr();



//-------------------------------------------------------------------//

Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> TPCmixture
(
    phaseChangeTwoPhaseMixture::New(U, phi)
);

//phaseChangeTwoPhaseMixture& TPCmixture = mixturePtr();

//volScalarField& alpha1(TPCmixture.alpha1());
//volScalarField& alpha2(TPCmixture.alpha2());

const dimensionedScalar& rho1 = TPCmixture->rho1();
const dimensionedScalar& rho2 = TPCmixture->rho2();


const dimensionedScalar& K1 = TPCmixture->K1();
const dimensionedScalar& K2 = TPCmixture->K2();

const dimensionedScalar& C1 = TPCmixture->C1();
const dimensionedScalar& C2 = TPCmixture->C2();

const volScalarField& TSatLocal = TPCmixture->TSatLocal();
const dimensionedScalar& pSat = TPCmixture->pSat();

//-------------------------------------------------------------------//

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + (scalar(1)-alpha1)*rho2,
    alpha1.boundaryField().types()
    //FluidThermo.rho()
);
rho.oldTime();
//-----------------------------------------------------------------//
Info<< "Generating field rhoPhi\n" << endl;

surfaceScalarField rhoPhi
(
    IOobject
    (
     	"rhoPhi",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*K1 + (scalar(1) - alpha1)*K2,
    alpha1.boundaryField().types()
);
//----------------------------------------------------------------//


 smoothInterfaceProperties interface(alpha1, U, TPCmixture());


// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, TPCmixture())
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;


setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

//-------------Flow work  -22.01.03-----------------//

Info<< "    Adding to flow work dpdt term\n" << endl;
 
 volScalarField dpdt
 (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            p.dimensions()/dimTime,
            0
        )
 );

//-------------------------------------------------//

#include "readTransportProperties.H"  //-----add 22.01.12

//-------------------mass transfer rate terms 22.01.03-----------------------------//
volScalarField mDotT
       (
       		IOobject
               (
                   "mDotT",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh,
               dimensionedScalar ("mDotT",dimDensity/dimTime,0)
       );

volScalarField mDotP
        (
        		IOobject
                (
                    "mDotP",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar ("mDotP",dimDensity/dimTime,0)
        );
volScalarField mDotAlpha
        (
            IOobject
                (
                    "mDotAlpha",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar ("mDotAlpha",dimDensity/dimTime,0)
         );

//---------------------------------------------------------------------------------// 
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createFvOptions.H"

