{
/*
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }
*/

    rAU = 1.0/UEqn.A();		//correcPhi 조건 제거
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));		//interface's RAU
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));	//openFOAM8에는 constraintHbyA가 새로 생성
									//-----22.01.13 추후 확인하기
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf)		//-----22.01.13 유도 방법 다름 확인하기
    );




//----------------------22.01.13 <2.2.0 버젼에는 없는 파트>---------------------//
    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }
//------------------------------------------------------------------------------//



    surfaceScalarField phig                         //source term?? by using H matrix?
    (
        (
            //mixture.surfaceTensionForce()               //correction needed
	    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)	    
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    Pair<tmp<volScalarField>> vDotP = TPCmixture->vDotP();

    const volScalarField& vDotcP = vDotP[0]();
    const volScalarField& vDotvP = vDotP[1]();







    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)

//--------------------22.01.13 압축성 기체의 압력의 의한 영향 고려 부분 ------------------//
//          - (vDotvP - vDotcP)*(TPCmixture->pSat() - rho*gh)
//          + fvm::Sp(vDotvP - vDotcP, p_rgh)
//----------------------------------------------------------------------------------------//

        );

        //p_rghEqn.setReference(pRefCell, pRefValue);
	p_rghEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));
//        p_rghEqn.solve();

//-------------------------------22.01.13------------------------------------//

	if (phaseChange){
		solve
			(
			 p_rghEqn == fvm::Sp(vDotcP-vDotvP,p_rgh)-(vDotcP-vDotvP)*(pSat-rho*gh)
			 //fvMatrix<scalar>::fvSolver::solve(p_rgh.select(pimple.finalInnerIter()))
			 //,mesh.solver(p_rgh.select(pimple.finalInnerIter()))
			);
	}

	else
	{
		solve
			(
			 p_rghEqn
			 //,mesh.solver(p_rgh.select(pimple.finalInnerIter()))
			);
	}
	
//--------------------------------------------------------------------------//

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            U = HbyA + rAU()*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

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
}
/*---------22.01.13
    if (!correctPhi)
    {
        rAU.clear();
    }
*/    
//------------------------------------22.01.13------------------------------------------//
//mass transfer rate 정의 for using at Energy Equation

    Pair<tmp<volScalarField> > mDot = TPCmixture->mDotP();
       const volScalarField& mDotcP = mDot[0]();
       const volScalarField& mDotvP = mDot[1]();

       mDotP = (mDotcP-mDotvP)*(p-pSat);

//--------------------------------------------------------------------------------------//
