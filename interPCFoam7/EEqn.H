Info<<"\n----------<Solve fluid hEqn>----------"<< endl;

{
    volScalarField& he = FluidThermo.he();
//    volScalarField& he = TPCmixture.he();
    const volScalarField& alphaDiffu = FluidThermo.alpha();

//---------------------------------22.01.13---------------------------------//

Pair<tmp<volScalarField> > vDotT = TPCmixture->vDotT();
									
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();	

	volScalarField K= K1*alpha1 + (1.0 - alpha1)*K2;
	volScalarField rhoC = rho1*C1*alpha1 + (1.0 - alpha1)*rho2*C2;


    volScalarField heT = he/rhoC;
//--------------------------------------------------------------------------//


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

//--------------------------------22.01.13-----------------------------------//

    if (phaseChange){
		solve
		 (
			EEqn
			==
			//fvm::Sp(vDotvT-vDotcT,T) - (vDotvT-vDotcT)*TSat
			//fvm::Sp(vDotcT-vDotvT,T)	//----22.01.13 T가 없어서 불가능
			fvm::Sp(vDotcT-vDotvT, he) 
			- (vDotcT-vDotvT)*TSatLocal	//
		 );
	 }
	 else
	 {
		 EEqn.solve();
	 }

//--------------------------------------------------------------------------//
    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);
    }
}


    FluidThermo.correct();
    //rad.correct();

    Info<< "Min/max T:" << min(FluidThermo.T()).value() << ' '
        << max(FluidThermo.T()).value() << endl;
        
//------------------------------22.01.13---------------------------------------//

Pair<tmp<volScalarField> > mDot = TPCmixture->mDotT();
	const volScalarField& mDotcT = mDot[0]();
	const volScalarField& mDotvT = mDot[1]();
		mDotT = (mDotcT-mDotvT)*(he-TSatLocal);

//-----------------------------------------------------------------------------//
Info<<"\n-----------------------------------------"<< endl;


