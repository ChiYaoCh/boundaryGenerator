  
    if (CoriolisOn)
    {
        fCoriolis = Omega ^ (Ug-U);
        forAll(fCoriolis,cellI)
        {
            fCoriolis[cellI].z() = 0.0;
        }
    }

    // Solve the momentum equation  
    tmp<fvVectorMatrix> UEqn
    (
      fvm::div(phi, U)     // convection
      + turbulence->divDevReff(U)  // turbulence modelling 
      - fCoriolis
      ==
	fvOptions(U)
    );

    UEqn().relax();
    fvOptions.constrain(UEqn());


   volVectorField thermoForce ("TherF",  fvc::reconstruct( (- ghf*fvc::snGrad((rhok))- fvc::snGrad(pd) )*mesh.magSf())); 

    if (simple.momentumPredictor())
    {
        solve
        (
            UEqn()
            ==
	      thermoForce

        );
        fvOptions.correct(U);
    }

    if(runTime.outputTime())
    {
     thermoForce.write();
     }

	
