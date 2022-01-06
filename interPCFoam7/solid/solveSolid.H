
Info<<"\n----------<Solve solid eEqn>----------"<< endl;

while (pimple2.loop())
{

volScalarField& e = wallThermo.he();

{
    while (pimple2.correctNonOrthogonal())
    {
    
     //solidT = wallThermo.T();	//wallThermo.T() value makes input of solidT field

         fvScalarMatrix eEqn
        (
	    fvm::ddt(wallThermo.rho(), e)//[J/s*m3] volumetric heat generation	//e[J/kg] from basicThermo
            								                                        //laplacian [1/m2]
	     -fvm::laplacian(wallThermo.alpha(), e, "laplacian(alpha,e)") //alpha [kg/m*s]
      	    ==
            fvOptions(rho, e)
        );

        eEqn.relax();

        fvOptions.constrain(eEqn);

        eEqn.solve();

        fvOptions.correct(e);
    }
}

wallThermo.correct();

Info<< "Min/max T:" << min(wallThermo.T()).value() << ' '
    << max(wallThermo.T()).value() << endl;

}

Info<<"\n-----------------------------------------"<< endl;
