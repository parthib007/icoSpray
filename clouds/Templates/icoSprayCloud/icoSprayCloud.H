/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::icoSprayCloud

Description
    Templated base class for spray cloud

    - sub-models:
      - atomization model
      - break-up model

\*---------------------------------------------------------------------------*/

#ifndef icoSprayCloud_H
#define icoSprayCloud_H

#include "icosprayCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CloudType>
class AtomizationModel;

template<class CloudType>
class BreakupModel;

/*---------------------------------------------------------------------------*\
                      Class SprayCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class icoSprayCloud
:
    public CloudType,
    public icosprayCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef icoSprayCloud<CloudType> icosprayCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<icoSprayCloud<CloudType>> cloudCopyPtr_;

        //- Average parcel mass
        scalar averageParcelMass_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        icoSprayCloud(const icoSprayCloud&);

        //- Disallow default bitwise assignment
        void operator=(const icoSprayCloud&);


protected:

    // Protected data

        // References to the cloud sub-models

            //- Atomization model
            autoPtr<AtomizationModel<icoSprayCloud<CloudType>>>
                atomizationModel_;

            //- Break-up model
            autoPtr<BreakupModel<icoSprayCloud<CloudType>>> breakupModel_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();

        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(icoSprayCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        icoSprayCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        icoSprayCloud(icoSprayCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        icoSprayCloud
        (
            const fvMesh& mesh,
            const word& name,
            const icoSprayCloud<CloudType>& c
        );


        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new icoSprayCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new icoSprayCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~icoSprayCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const icoSprayCloud& cloudCopy() const;

            //- Return const-access to the average parcel mass
            inline scalar averageParcelMass() const;


        // Check

            //- Penetration for fraction [0-1] of the current total mass
            inline scalar penetration(const scalar fraction) const;


            // Sub-models

                //- Return const-access to the atomization model
                inline const AtomizationModel<icoSprayCloud<CloudType>>&
                    atomization() const;

                //- Return reference to the atomization model
                inline AtomizationModel<icoSprayCloud<CloudType>>& atomization();

                //- Return const-access to the breakup model
                inline const BreakupModel<icoSprayCloud<CloudType>>&
                    breakup() const;

                //- Return reference to the breakup model
                inline BreakupModel<icoSprayCloud<CloudType>>& breakup();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the spray (inject, move)
            void evolve();

            //- Inject parcels
            void inject
            (
    		const fvMesh& mesh,
    		const scalar dt,
		vector position,
		const scalar diameter,
		const vector velocity
            );

        // I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "icoSprayCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "icoSprayCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
