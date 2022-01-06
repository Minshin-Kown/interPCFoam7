scalar DiNum = -great;

{
    // Note: do not use setRegionSolidFields.H to avoid double registering Cp
    //#include "setRegionSolidFields.H"
   // const solidThermo& wallThermo = solidThermo;

    tmp<volScalarField> magKappa;
	
    magKappa = wallThermo.kappa();

    tmp<volScalarField> tcp = wallThermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> trho = wallThermo.rho();
    const volScalarField& rho = trho();

    DiNum = max
    (
        solidRegionDiffNo
        (
            mesh2,
            runTime,
            rho*cp,
            magKappa()
        ),
        DiNum
    );
}
