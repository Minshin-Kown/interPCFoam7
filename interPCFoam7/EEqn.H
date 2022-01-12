Info<<"\n----------<Solve fluid hEqn>----------"<< endl;

{
    volScalarField& he = FluidThermo.he();
    const volScalarField& alphaDiffu = FluidThermo.alpha();

{
    while (pimple.correctNonOrthogonal())
    {
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(rhoPhi, he) 		//phi is including rho [kg/s]
      //+ fvc::ddt(rho, KE) + fvc::div(rhoPhi, KE)			//K is (1/2)u^2; K.E. term
      + (
         //   he.name() == "e"
         // ?
	 // fvc::div						//internal energy term
         //   (
         //       fvc::absolute(rhoPhi/fvc::interpolate(rho), U), 	
         //       p,
         //       "div(phiv,p)"
	 // )
          -dpdt						//flow work term (Pdv devided by rho)
        )
      //+ thermophysicalTransport.divq(he)
         -fvm::laplacian (fvc::interpolate(alphaDiffu), he, "(alphaDiffu,he)")

     ==
        rho*(U&g)
      //+ rad.Sh(thermo, he)					//radiation term = off
      //+ reaction.Qdot()					//reaction term = off
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);
    }
}


    FluidThermo.correct();
    //rad.correct();

    Info<< "Min/max T:" << min(FluidThermo.T()).value() << ' '
        << max(FluidThermo.T()).value() << endl;
Info<<"\n-----------------------------------------"<< endl;


