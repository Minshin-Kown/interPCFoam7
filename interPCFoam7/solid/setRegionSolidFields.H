

//solidNoLoopControl& pimple = pimples.solid();

/*
tmp<volScalarField> trho = wallThermo.rho();
const volScalarField& solidRho = trho();

tmp<volScalarField> tCv = wallThermo.Cv();
const volScalarField& solidCv = tCv();


tmp<volSymmTensorField> taniAlpha;
if (!wallThermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCv = wallThermo.Kappa()/Cv;

    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.R().transformVector(tkappaByCv());
    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );

}


volScalarField& e = wallThermo.he();

const volScalarField& betav = betavSolid();


fv::options& fvOptions = solidHeatSources[i];

solidNoLoopControl& pimple = pimples.solid(i);
*/

pimpleControl pimple(mesh2);
