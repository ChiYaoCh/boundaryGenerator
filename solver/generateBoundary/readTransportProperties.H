    Info<< "\nReading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );


    dimensionedScalar nu( transportProperties.lookup("nu"));

    // Thermal expansion coefficient [1/K]
//    dimensionedScalar beta(transportProperties.lookup("beta"));

    // Reference temperature [K]
    dimensionedScalar TRef(transportProperties.lookup("TRef"));

    // Laminar Prandtl number
    dimensionedScalar Pr(transportProperties.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(transportProperties.lookup("Prt"));

    // Specific Heat capacity
    dimensionedScalar Cp(transportProperties.lookup("Cp"));

    // gravitation acceleration
    dimensionedVector g(transportProperties.lookup("g"));

    // Coriolis coefficient
    dimensionedVector Omega(transportProperties.lookup("Omega"));

    //  Geostrophic Velocity
    dimensionedVector Ug(transportProperties.lookup("Ug") );

    //  Read centerline velocity for channel simulations
    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );

    dimensionedScalar magUbar = mag(Ubar);
    vector flowDirection = (Ubar/magUbar).value();

