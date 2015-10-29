if (tempEqnOn)
{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::Sp(fvc::div(phi), T)
      - fvm::laplacian(alphaEff, T)
    );

    TEqn.relax();
    TEqn.solve();

//    rhok = 1.0- (T - TRef)/TRef;
    rhok = (TRef - T)/TRef;       	
}

