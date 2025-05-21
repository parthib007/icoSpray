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

#include "icoSprayParcel.H"
#include "AtomizationModel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::icoSprayParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    ParcelType::setCellValues(td, dt, celli);
}


template<class ParcelType>
template<class TrackData>
void Foam::icoSprayParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    ParcelType::cellValueSourceCorrection(td, dt, celli);
}


template<class ParcelType>
template<class TrackData>
void Foam::icoSprayParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    //typedef typename TrackData::cloudType::icosprayCloudType icosprayCloudType;

    // Check if parcel belongs to liquid core
    if (liquidCore() > 0.5)
    {
        // Liquid core parcels should not experience coupled forces
           td.cloud().forces().setCalcCoupled(false);
    }

    ParcelType::calc(td, dt, celli);

    if (td.keepParticle)
    {
        scalar d1 = this->d();

        if (liquidCore() > 0.5)
        {
            calcAtomization(td, dt, celli);

            // Preserve the total mass/volume by increasing the number of
            // particles in parcels due to breakup
            scalar d2 = this->d();
            this->nParticle() *= pow3(d1/d2);
            //calcBreakup(td, dt, celli); //if blob breakup model
        }
        else
        {
            calcBreakup(td, dt, celli);
        }
    }

    // Restore coupled forces
    td.cloud().forces().setCalcCoupled(true);
}


template<class ParcelType>
template<class TrackData>
void Foam::icoSprayParcel<ParcelType>::calcAtomization
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    typedef typename TrackData::cloudType::icosprayCloudType icosprayCloudType;
    const AtomizationModel<icosprayCloudType>& atomization =
        td.cloud().atomization();
    //Info<<"\nPrimaryBreakupBegins"<<endl;
    // Store average contious phase density
    scalar rhoAv = this->rhoc();

    scalar soi = td.cloud().injectors().timeStart();
    scalar currentTime = td.cloud().db().time().value();
    const vector& pos = this->position();
    const vector& injectionPos = this->position0();

    // Disregard the continous phase when calculating the relative velocity
    // (in line with the deactivated coupled assumption)
    scalar Urel = mag(this->U());

    scalar t0 = max(0.0, currentTime - this->age() - soi);
    scalar t1 = min(t0 + dt, td.cloud().injectors().timeEnd() - soi);

    // This should be the vol flow rate from when the parcel was injected
    scalar volFlowRate = td.cloud().injectors().volumeToInject(t0, t1)/dt;

    atomization.update
    (
        dt,
        this->d(),
        this->liquidCore(),
        this->tc(),
        this->rho(),
        mu_,
        sigma_,
        volFlowRate,
        rhoAv,
        Urel,
        pos,
        injectionPos,
        td.cloud().pAmbient(),
        td.cloud().rndGen()
    );
}


template<class ParcelType>
template<class TrackData>
void Foam::icoSprayParcel<ParcelType>::calcBreakup
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::forceType forceType;

    const parcelType& p = static_cast<const parcelType&>(*this);
    const forceType& forces = td.cloud().forces();

    if (td.cloud().breakup().solveOscillationEq())
    {
        solveTABEq(td, dt);
    }

    // Calculate average continous phase properties
    scalar rhoAv = this->rhoc();
    scalar muAv = this->muc();
    vector Urel = this->U()- this->Uc();
    scalar Urmag = mag(Urel);
    scalar Re = this->Re(this->U(), this->d(), rhoAv, muAv);

    const scalar mass = p.mass();
    const forceSuSp Fcp = forces.calcCoupled(p, dt, mass, Re, muAv);
    const forceSuSp Fncp = forces.calcNonCoupled(p, dt, mass, Re, muAv);

    if(Fcp.Sp() != 0 || Fncp.Sp() != 0)
    {
    this->tMom() = mass/(Fcp.Sp() + Fncp.Sp());  //check this
    }

    const vector g = td.cloud().g().value();

    scalar parcelMassChild = 0.0;
    scalar dChild = 0.0;

    if
    (
        td.cloud().breakup().update
        (
            dt,
            g,
            this->d(),
            this->tc(),
            this->ms(),
            this->nParticle(),
            this->KHindex(),
            this->y(),
            this->yDot(),
            this->d0(),
            this->rho(),
            mu_,
            sigma_,
            this->U(),
            rhoAv,
            muAv,
            Urel,
            Urmag,
            this->tMom(),
            dChild,
            parcelMassChild
        )
    )
    {
        scalar Re = rhoAv*Urmag*dChild/muAv;

        // Add child parcel as copy of parent
        icoSprayParcel<ParcelType>* child = new icoSprayParcel<ParcelType>(*this);
        child->origId() = this->getNewParticleID();
        child->d() = dChild;
        child->d0() = dChild;
        const scalar massChild = child->mass();
        //child->mass0() = massChild;
        child->nParticle() = parcelMassChild/massChild;

        const forceSuSp Fcp =
            forces.calcCoupled(*child, dt, massChild, Re, muAv);
        const forceSuSp Fncp =
            forces.calcNonCoupled(*child, dt, massChild, Re, muAv);

        child->age() = 0.0;
        child->liquidCore() = 0.0;
        child->KHindex() = 1.0;
        child->y() = td.cloud().breakup().y0();
        child->yDot() = td.cloud().breakup().yDot0();
        child->tc() = 0.0;
        child->ms() = -GREAT;
        child->injector() = this->injector();
    	if(Fcp.Sp() != 0 || Fncp.Sp() != 0)
    	{
        child->tMom() = massChild/(Fcp.Sp() + Fncp.Sp());
	}	
        child->user() = 0.0;
        child->setCellValues(td, dt, celli);

        td.cloud().addParticle(child);
    }
}

template<class ParcelType>
template<class TrackData>
void Foam::icoSprayParcel<ParcelType>::solveTABEq
(
    TrackData& td,
    const scalar dt
)
{
    const scalar& TABCmu = td.cloud().breakup().TABCmu();
    const scalar& TABtwoWeCrit = td.cloud().breakup().TABtwoWeCrit();
    const scalar& TABComega = td.cloud().breakup().TABComega();

    scalar r = 0.5*this->d();
    scalar r2 = r*r;
    scalar r3 = r*r2;
    
    // Inverse of characteristic viscous damping time
    scalar rtd = 0.5*TABCmu*mu_/(this->rho()*r2);

    // Oscillation frequency (squared)
    scalar omega2 = TABComega*sigma_/(this->rho()*r3) - rtd*rtd;
    //Info << "I am now here!" << endl;
    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar rhoc = this->rhoc();
        scalar We = this->We(this->U(), r, rhoc, sigma_)/TABtwoWeCrit;

        // Initial values for y and yDot
        scalar y0 = this->y() - We;
        scalar yDot0 = this->yDot() + y0*rtd;

        // Update distortion parameters
        scalar c = cos(omega*dt);
        scalar s = sin(omega*dt);
        scalar e = exp(-rtd*dt);

        this->y() = We + e*(y0*c + (yDot0/omega)*s);
        this->yDot() = (We - this->y())*rtd + e*(yDot0*c - omega*y0*s);
    }
    else
    {
        // Reset distortion parameters
        this->y() = 0;
        this->yDot() = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::icoSprayParcel<ParcelType>::icoSprayParcel(const icoSprayParcel<ParcelType>& p)
:
    ParcelType(p),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_)
{}


template<class ParcelType>
Foam::icoSprayParcel<ParcelType>::icoSprayParcel
(
    const icoSprayParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "icoSprayParcelIO.C"


// ************************************************************************* //
