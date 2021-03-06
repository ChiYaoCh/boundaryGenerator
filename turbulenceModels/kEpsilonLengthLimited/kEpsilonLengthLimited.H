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


Class
    Foam::incompressible::RASModels::kEpsilonLengthLimited

Description
    k-epsilon turbulence model for incompressible flows with
    Apsley and Castro correction for finite upper limit of mixing length.

    The default model coefficients correspond to the following:
    \verbatim
        kEpsilonLengthLimitedCoeffs
        {
            Cmu         0.09;
            C1          = C1 + (C2-C1)*lm/lmax
            C2          1.92;
            sigmaEps    1.3;
            lmax        1.0  maximum limited mixing length
        }
    \endverbatim

SourceFiles
    kEpsilonLengthLimited.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonLengthLimited_H
#define kEpsilonLengthLimited_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonLengthLimited Declaration
\*---------------------------------------------------------------------------*/

class kEpsilonLengthLimited
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C4_;
            dimensionedScalar sigmaEps_;

        // Fields
            word UgName_;
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;


public:

    //- Runtime type information
    TypeName("kEpsilonLengthLimited");

    // Constructors

        //- Construct from components
        kEpsilonLengthLimited
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~kEpsilonLengthLimited()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;


        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
