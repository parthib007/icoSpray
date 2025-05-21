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

\*---------------------------------------------------------------------------*/

#include "icoSprayCloud.H"
#include "AtomizationModel.H"
#include "BreakupModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::icoSprayCloud<CloudType>::setModels()
{
    atomizationModel_.reset
    (
        AtomizationModel<icoSprayCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    breakupModel_.reset
    (
        BreakupModel<icoSprayCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::icoSprayCloud<CloudType>::cloudReset
(
    icoSprayCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    atomizationModel_.reset(c.atomizationModel_.ptr());
    breakupModel_.reset(c.breakupModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::icoSprayCloud<CloudType>::icoSprayCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
    icosprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(0.0),
    atomizationModel_(NULL),
    breakupModel_(NULL)
{
    if (this->solution().active())
    {
        setModels();

        averageParcelMass_ = this->injectors().averageParcelMass();

        if (readFields)
        {
            parcelType::readFields(*this);
        }

        Info << "Average parcel mass: " << averageParcelMass_ << endl;
    }

    if (this->solution().resetSourcesOnStartup())
    {
        CloudType::resetSourceTerms();
    }
}


template<class CloudType>
Foam::icoSprayCloud<CloudType>::icoSprayCloud
(
    icoSprayCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    icosprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(c.averageParcelMass_),
    atomizationModel_(c.atomizationModel_->clone()),
    breakupModel_(c.breakupModel_->clone())
{}


template<class CloudType>
Foam::icoSprayCloud<CloudType>::icoSprayCloud
(
    const fvMesh& mesh,
    const word& name,
    const icoSprayCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    icosprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(0.0),
    atomizationModel_(NULL),
    breakupModel_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::icoSprayCloud<CloudType>::~icoSprayCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::icoSprayCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);
    const scalar rhop=readScalar(this->particleProperties().subDict("constantProperties").lookup("rhop"));
    const scalar mup=readScalar(this->particleProperties().subDict("constantProperties").lookup("mup"));
    const scalar sigmap=readScalar(this->particleProperties().subDict("constantProperties").lookup("sigmap"));
    // override from constantProperties 
    parcel.rho() = rhop;
    parcel.mu() = mup;
    parcel.sigma() = sigmap;
}


template<class CloudType>
void Foam::icoSprayCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    // store the injection position and initial drop size
    parcel.position0() = parcel.position();
    parcel.d0() = parcel.d();

    parcel.y() = breakup().y0();
    parcel.yDot() = breakup().yDot0();

    parcel.liquidCore() = atomization().initLiquidCore();
}


template<class CloudType>
void Foam::icoSprayCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<icoSprayCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::icoSprayCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::icoSprayCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<icoSprayCloud<CloudType>> td(*this);

        this->solve(td);
    }
}


template<class CloudType>
void Foam::icoSprayCloud<CloudType>::info()
{
    CloudType::info();
    scalar d32 = 1.0e+6*this->Dij(3, 2);
    scalar d10 = 1.0e+6*this->Dij(1, 0);
    scalar dMax = 1.0e+6*this->Dmax();
    scalar pen = this->penetration(0.95);

    Info << "    D10, D32, Dmax (mu)             = " << d10 << ", " << d32
         << ", " << dMax << nl
         << "    Liquid penetration 95% mass (m) = " << pen << endl;
}

template<class CloudType>
void Foam::icoSprayCloud<CloudType>::inject
(
    const fvMesh& mesh,
    const scalar dt,
    vector position,
    const scalar diameter,
    const vector velocity
)
{
    label celli;
    label tetFacei;
    label tetPti;

    volVectorField cellCentres = mesh.C();

        typename parcelType::template
            TrackingData<icoSprayCloud<CloudType>> td(*this);

    // Find the celli, tetFacei and tetPti for point position
    mesh.findCellFacePt
    (
        position,
        celli,
        tetFacei,
        tetPti
    );

    label proci = -1;

    if (celli >= 0)
    {
        proci = Pstream::myProcNo();
    }

    reduce(proci, maxOp<label>());
	Info << proci << endl;

    // Ensure that only one processor attempts to insert this Parcel
    if (proci != Pstream::myProcNo())
    {
        celli = -1;
        tetFacei = -1;
        tetPti = -1;
    }

    // Last chance - find nearest cell and try that one - the point is
    // probably on an edge
    if (proci == -1)
    {
        celli = mesh.findNearestCell(position);

        if (celli >= 0)
        {
            position += 1e-3*(cellCentres[celli] - position);

            mesh.findCellFacePt
            (
                position,
                celli,
                tetFacei,
                tetPti
            );

            if (celli > 0)
            {
                proci = Pstream::myProcNo();
            }
        }

        reduce(proci, maxOp<label>());

        if (proci != Pstream::myProcNo())
        {
            celli = -1;
            tetFacei = -1;
            tetPti = -1;
        }
    }

    if (celli > -1)
    {
	parcelType* pPtr =
                        new parcelType(mesh, position, celli, tetFacei, tetPti);
        pPtr->origId() = pPtr->getNewParticleID();
	pPtr->d0()=diameter;
	pPtr->U()=velocity;
	pPtr->position0()=position;
	pPtr->d()=diameter;
	pPtr->y()=0;
	pPtr->yDot()=0;
	pPtr->tc()=0;
        pPtr->age() = 0.0;
        pPtr->liquidCore() = 0.0;
        pPtr->KHindex() = 1.0;
        pPtr->ms() = -GREAT;
        pPtr->nParticle() = 1.0;
        pPtr->user() = 0.0;
    	pPtr->rho() = readScalar(this->particleProperties().subDict("constantProperties").lookup("rhop"));;
    	pPtr->mu() = readScalar(this->particleProperties().subDict("constantProperties").lookup("mup"));
    	pPtr->sigma() = readScalar(this->particleProperties().subDict("constantProperties").lookup("sigmap"));
        pPtr->setCellValues(td, dt, celli);
        td.cloud().addParticle(pPtr);
	//td.cloud().averageParcelMass() == 0.0;
    }
}

// ************************************************************************* //
