/*---------------------------------------------------------------------------*\
         _      _                 |
        | |    | |                | OLDEV: CFD tools developed by
   ___  | |  __| |  ___ __   __   |        Fraunhofer IWES and ForWind
  / _ \ | | / _` | / _ \\ \ / /   |        in Oldenburg, Germany
 | (_) || || (_| ||  __/ \ V /    |
  \___/ |_| \__,_| \___|  \_/     | http://www.iwes.fraunhofer.de
                                  | http://www.forwind.de
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OLDEV and it is based on OpenFOAM.

    OLDEV is owned by the Fraunhofer IWES and ForWind Oldenburg CFD groups.
    It is not allowed to use it outside the scope of Fraunhofer IWES or
    ForWind Oldenburg CFD projects. It is forbidden to use, modify, copy or
    redistribute OLDEV files for any other purpose.

    OLDEV is developed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.

    OpenFOAM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.



\*---------------------------------------------------------------------------*/

#include "kEpsilonLengthLimited.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilonLengthLimited, 0);
addToRunTimeSelectionTable(RASModel, kEpsilonLengthLimited, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilonLengthLimited::kEpsilonLengthLimited
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
        
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            coeffDict_,
            0.33
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),

   UgName_
   (
    coeffDict_.lookupOrDefault<word>("UgName","Ug")
   ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kEpsilonLengthLimited::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kEpsilonLengthLimited::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kEpsilonLengthLimited::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

tmp<fvVectorMatrix> kEpsilonLengthLimited::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool kEpsilonLengthLimited::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsilonLengthLimited::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }


    const dimensionedVector& Ug_ = U_.db().lookupObject<IOdictionary >("transportProperties").lookup("Ug");

    const dimensionedVector& omega_ = U_.db().lookupObject< IOdictionary >("transportProperties").lookup("Omega");

    volScalarField G(GName(), nut_*2*magSqr(symm(fvc::grad(U_))));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    dimensionedScalar smallS_("smallS", dimless/dimTime, 1e-10);
    dimensionedScalar smallU_("smallU", dimLength/dimTime, 1e-10);
//    volScalarField lmax_("lmax", 0.00027*(mag(Ug_)+smallU_)/(mag(omega_)+smallS_));
    dimensionedScalar lmax_("lmax", 0.00027*(mag(Ug_)+smallU_)/(mag(omega_)+smallS_));
    volScalarField Lboun = pow(Cmu_,0.75)*pow(k_,1.5)/epsilon_ ;
    volScalarField C1prime_ = C1_ + (C2_ - C1_)*Lboun/lmax_;

// change here C1 following suggestion of Apsley and Castro 97:
    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::Sp(fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
         C1prime_*G*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)

    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
