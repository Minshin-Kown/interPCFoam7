
//------------------------------22.01.07----------------------//
surfaceScalarField phic(mag(phi/mesh.magSf()));
phic = min(interface.cAlpha()*phic, max(phic)); 
//-------------------------------------------------------------//

//----------------------------22.01.12-------------------------//
dimensionedScalar totalDeltaT = runTime.deltaT();
//-------------------------------------------------------------//

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), 0)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
//        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;		//----22.01.12
    }

//    rhoPhi = rhoPhiSum;	//----22.01.12
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + (scalar(1)-alpha1)*rho2;
