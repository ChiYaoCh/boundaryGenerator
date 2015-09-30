
    Info<< "Creating/Reading deviation from hydrostatic pressure field, pd..." << endl;

    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    bool tempEqnOn = mesh.solutionDict().subDict("options").lookupOrDefault<bool>("tempEqnOn", false);
    bool CoriolisOn = mesh.solutionDict().subDict("options").lookupOrDefault<bool>("CoriolisOn", false);

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    singlePhaseTransportModel laminarTransport(U, phi);
    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );


    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar ("rhok", dimless, 0.0)
    );
    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf()); 

    Info<< "Reading thermophysical properties\n" << endl;

    /* Create and read in initial "potential temperature" field	*/
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField T_abs
    (
        IOobject
        (
            "T_abs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T + gh/Cp
    );


    // kinematic turbulent thermal thermal conductivity m2/s
    Info<< "Reading field alphat\n" << endl;
    volScalarField alphat
    (
        IOobject
        (
            "alphat",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating the Coriolis force vector, fCoriolis..." << endl;
    volVectorField fCoriolis
    (
        IOobject
        (
            "fCoriolis",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("fCoriolis",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
    );
