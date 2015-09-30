    volScalarField rAU("rAU", 1.0/UEqn().A());

{
    surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

    U = rAU*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(U) & mesh.Sf();

    fvOptions.makeRelative(phi);
    adjustPhi(phi, U, pd);

    surfaceScalarField buoyancyPhi(rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

    phi -= buoyancyPhi;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        pd.boundaryField(),
        (
            phi.boundaryField()
          - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );


    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAUf, pd) == fvc::div(phi)
        );

   //     pdEqn.setReference(pdRefCell, getRefCellValue(pd, pdRefCell));

        pdEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi -= pdEqn.flux();

            // Explicitly relax pressure for momentum corrector
            pd.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rAU*fvc::reconstruct((buoyancyPhi + pdEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

}
