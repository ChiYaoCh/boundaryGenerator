    dimensionedScalar gradPd
    (
        "gradPd",
        dimensionSet(0, 1, -2, 0, 0),
        0.0
    );

    IFstream gradPdFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"gradPd.raw"
    );

    if(gradPdFile.good())
    {
        gradPdFile >> gradPd;
        Info<< "Reading average pressure (pd) gradient..." << endl << endl;
    }
    else
    {
        Info<< "Initializing with 0 pressure (pd) gradient..." << endl << endl;
    };
