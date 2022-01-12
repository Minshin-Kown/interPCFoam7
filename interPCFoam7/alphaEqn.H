{
    #include "alphaScheme.H"

//------------------------22.01.10-----------------------------------------//

//add alpha"R"Scheme
	
static const word divAlpharName("div(phirb,alpha)");

const word alpharScheme(mesh.divScheme(divAlpharName)[1].wordToken());

ITstream compressionScheme2
(
    compressionSchemes.found(alpharScheme)
  ? mesh.divScheme(divAlpharName)
  : ITstream
    (
        divAlpharName,
        tokenList
        {
            word("Gauss"),
            word("interfaceCompression"),
            alpharScheme,
            alphaControls.lookup<scalar>("cAlpha")
        }
    )
);
//---------------------------------------------------------------------------//

    Pair<tmp<volScalarField>> vDotAlphal = TPCmixture->vDotAlphal();
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
   // const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    //tmp<surfaceScalarField> talphaPhi;
/*---------------------------------------22.01.12--------------------------//
    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        talphaPhi = alpha1Eqn.flux();
    }
----------------------------------------------------------------------------------*/
    volScalarField alpha10("alpha10", alpha1);
//--------------------22.01.07 chtFlowboiling 시작부분-------------------------------------------//

//------------------------------22.01.07---------------------------------------------------------//
	surfaceScalarField phir("phir", phic*interface.nHatf());
//-----------------------------------------------------------------------------------------------//

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
 	surfaceScalarField  phiAdd = -fvc::flux(-phir, scalar(1)-alpha1, divAlpharName); //-----22.01.12 tmp 안됨

        surfaceScalarField talphaPhiCorr
	//surfaceScalarField phiAlpha	//오류 발생
	(
            fvc::flux
            (
                phi,
                alpha1,
                compressionScheme.rewind()
            )
//------------------------------22.01.12--------------------------------------------------------//
	    + fvc::flux
	    (
         	phiAdd,		//overloading part
		alpha1,
		compressionScheme2.rewind()
	    )
        );

       volScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            vDotvAlphal - vDotcAlphal
        );


        volScalarField Su
        (   
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*alpha1
          + vDotcAlphal
        );

	if(phaseChange)
	{
		MULES::explicitSolve
			(
			 geometricOneField(),
			 alpha1,
			 phi,
			 //phiAlpha(),
			 talphaPhiCorr,
			 Sp,
			 Su,
			 oneField(),
			 zeroField()
			);
	}
	else
	{
		MULES::explicitSolve
			(
			 geometricOneField(),
			 alpha1, 
			 phi, 
			 talphaPhiCorr,
			 //phiAlpha, 
			 oneField(), 
			 zeroField()
			);
	}

//----------------------------------------------------------------------------------------------//

/*---------------------------------22.01.12---------------------------------------------------//
        if (MULESCorr)
        {
            talphaPhiCorr.ref() -= talphaPhi();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhi(),
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (
                    divU*(alpha10 - alpha100)
                  - vDotvmcAlphal*alpha10
                )(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talphaPhi.ref() += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                talphaPhi.ref() += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (divU*alpha1 + vDotcAlphal)(),
                oneField(),
                zeroField()
            );

            talphaPhi = talphaPhiCorr;
        }
----------------------------------------------------------------------------------------------*/
       // alpha2 = 1.0 - alpha1;
       rhoPhi +=
	       (runTime.deltaT()/totalDeltaT)
	       *(talphaPhiCorr*(rho1 - rho2) + phi*rho2);
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
