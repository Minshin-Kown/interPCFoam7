Info<< "--------------------------------\n" << endl;
Info<< "\nCreate solid fields\n" << endl;
/*
Info<< "Adding field solid T\n" << endl;
volScalarField solidT			//2 for solid
(
    IOobject
    (
        "T",				//read name from library
        runTime.timeName(),
        mesh2,				//2 for solid
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh2				//2 for solid
);
*/
Info<< "Adding field solid p\n" << endl;
volScalarField solidP			
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh2,	
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh2				
);

Info<< "Reading solid thermophysical properties\n" << endl;

autoPtr<solidThermo> mySolidThermoPtr (solidThermo::New(mesh2));

solidThermo& wallThermo = mySolidThermoPtr();



/*
Info<< "Adding field solid alpha\n" << endl;  //for porous material
volScalarField solidAlpha
(
     IOobject
     (
         "alpha",
         runTime.timeName(),
         mesh2,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
     wallThermo.alpha()
);
 */


Info<< "Adding field solid rho\n" << endl;
volScalarField solidRho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh2,
        IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    wallThermo.rho()

);

/*
Info<< "Adding field KbyCv\n" << endl;

volScalarField tkappaByCp
(
     IOobject
         (
         "KbyCp",
         runTime.timeName(),
         mesh2,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
         ),
         wallThermo.kappa()/wallThermo.Cp()
);
*/

/*
volScalarField betavSolid
(
       IOobject
       (
            "betavSolid",
            runTime.timeName(),
   	    mesh2,
       	    IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       mesh2,dimensionedScalar(dimless, scalar(1))
);
*/

//const volScalarField& betav = betavSolid;

