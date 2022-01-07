/*
surfaceScalarField rhoPhi
     (
        IOobject
           (
             "rhoPhi",
             runTime.timeName(),
             mesh
           ),
         mesh,
         dimensionedScalar(dimMass/dimTime, 0)
     );
*/

const dictionary& alphaControls = mesh.solverDict(alpha1.name());	//fvSolution에서 alpha1.name 읽어와 별칭 alphaControls

const label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));	//alphaControls의 nAlphaCoor

const label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));

const bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
// const bool alphaApplyPrevCorr
// (
//     alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false)
// );
