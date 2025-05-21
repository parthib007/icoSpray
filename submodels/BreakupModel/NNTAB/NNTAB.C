/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANNNTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "NNTAB.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NNTAB<CloudType>::NNTAB
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict, owner, typeName),
    y0_(0),
    yDot0_(0),
    WeCrit1_(12),
    WeCrit2_(35),
    WeCrit3_(80),
    NNindex_(1.0),
    NP_(50.0),
    C_(0.5),
    Factor1_(1.2),
    Factor2_(1.0),
    SMDCalcMethod_(word("method1"))
{
	if(!this->defaultCoeffs(true))
	{
		this->coeffDict().lookup("y0");
		this->coeffDict().lookup("yDot0");
		this->coeffDict().lookup("WeCrit1");
		this->coeffDict().lookup("WeCrit2");
		this->coeffDict().lookup("WeCrit3");
		this->coeffDict().lookup("NNindex");
		this->coeffDict().lookup("NP");
		this->coeffDict().lookup("C");
		this->coeffDict().lookup("Factor1");
		this->coeffDict().lookup("Factor2");
		this->coeffDict().lookup("SMDCalculationMethod");
	}
    // calculate the inverse function of the Rossin-Rammler Distribution
    const scalar xx0 = 12.0;
    const scalar rrd100 =
        1.0/(1.0 - exp(-xx0)*(1.0 + xx0 + sqr(xx0)/2.0 + pow3(xx0)/6.0));

    forAll(rrd_, n)
    {
        scalar xx = 0.12*(n + 1);
        rrd_[n] =
            (1.0 - exp(-xx)*(1.0 + xx + sqr(xx)/2.0 + pow3(xx)/6.0))*rrd100;
    }

    if (SMDCalcMethod_ == "method1")
    {
        SMDMethod_ = method1;
    }
    else if (SMDCalcMethod_ == "method2")
    {
        SMDMethod_ = method2;
    }
    else
    {
        SMDMethod_ = method2;
        WarningInFunction
            << "Unknown SMDCalculationMethod. Valid options are "
            << "(method1 | method2). Using method2" << endl;
    }
}


template<class CloudType>
Foam::NNTAB<CloudType>::NNTAB(const NNTAB<CloudType>& bum)
:
    BreakupModel<CloudType>(bum),
    y0_(bum.y0_),
    yDot0_(bum.yDot0_),
    WeCrit1_(bum.WeCrit1_),
    WeCrit2_(bum.WeCrit2_),
    WeCrit3_(bum.WeCrit3_),
    NNindex_(bum.NNindex_),
    NP_(bum.NP_),
    C_(bum.C_),
    Factor1_(bum.Factor1_),
    Factor2_(bum.Factor2_),
    SMDCalcMethod_(bum.SMDCalcMethod_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NNTAB<CloudType>::~NNTAB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NNTAB<CloudType>::update
(
    const scalar dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar d0,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const vector& U,
    const scalar rhoc,
    const scalar muc,
    const vector& Urel,
    const scalar Urmag,
    const scalar tMom,
    scalar& dChild,
    scalar& massChild
)
{
    bool addParcel = false;

    cachedRandom& rndGen = this->owner().rndGen();

    scalar r = 0.5*d;
    scalar r2 = r*r;
    scalar r3 = r*r2;

    scalar semiMass = nParticle*pow3(d);

    // droplet deformation characteristic time
    scalar tChar = d/Urmag*sqrt(rho/rhoc);
    scalar tFirst = 0.5*tChar;
    scalar tSecond = 0;
    scalar tCharSecond = 0;

    // update the droplet characteristic time
    tc += dt;

    // Calculating theoritical critical weber number
    scalar weGas = rhoc*sqr(Urmag)*d/sigma;
    scalar weLiquid = rho*sqr(Urmag)*d/sigma;
    scalar reLiquid = rho*Urmag*d/mu;
    scalar ohnesorge = sqrt(weLiquid)/(reLiquid + VSMALL);
    scalar densitycorr = 1+(rhoc/(rho-rhoc));
    scalar ustar = 3*(muc/mu)*(1+log10(1000*rhoc/rho));
    scalar Zhi1 = 2.25*(1-ustar)*(1-ustar);
    scalar Zhi2 = 4.00*(1-ustar)*(1-ustar);
    scalar Zhi3 = 4.00*(1-ustar)*(1-ustar); 
    scalar wethcr1 = WeCrit1_*(1+1.4*pow(Zhi1,NNindex_-1)*ohnesorge)*pow(densitycorr,2);
    scalar wethcr2 = WeCrit2_*(1+1.4*pow(Zhi2,NNindex_-1)*ohnesorge)*pow(densitycorr,1);
    scalar wethcr3 = WeCrit3_*(1+1.4*pow(Zhi3,NNindex_-1)*ohnesorge)*pow(densitycorr,1); 
    scalar sr = 1;
    scalar Tbal = (tFirst*(1+ohnesorge))/tChar;

    if ( (weGas > wethcr1) && (weGas < wethcr2) ) //bag breakup
    {
	sr = 3*0.25*(1+log10(rho/rhoc));
        tCharSecond = 6*pow((weGas - 12), 0.25);
    	yDot = (pow(0.5*sr,2)-(16/weGas))*Tbal-16*(ohnesorge/sqrt(weGas));
    	y = 1.5;
    }
    else if ( (weGas > wethcr2) && (weGas < wethcr3) ) //multimode breakup 
    {
	sr = 3.25*0.25*(1+log10(rho/rhoc));
        tCharSecond = 2.45*pow((weGas - 12), -0.25);
    	yDot = (pow(0.5*sr,2)-(16/weGas))*Tbal-16*(ohnesorge/sqrt(weGas));
    	y = 2.0;
    }
    else if ( (weGas > wethcr3) ) //sheet thinning breakup 
    {
	sr = 3.5*0.25*(1+log10(rho/rhoc));
        tCharSecond = 14.1*pow((weGas - 12), 0.25);
    	yDot = (pow(0.5*sr,2)-(16/weGas))*Tbal-16*(ohnesorge/sqrt(weGas));
    	y = 2.0;
    }
    else //no breakup
    {
	sr = 2.77;
        y = 0;
        yDot = 0;
    }

    tSecond = tCharSecond*tChar;
    scalar tBreakUP = tFirst + tSecond;
    tBreakUP = tBreakUP*(1+ohnesorge);

    // update droplet size
    if ( (tc > tBreakUP) && (weGas > wethcr1) )
    {
	scalar rs = r/(1.0 + (4.0/3.0)*sqr(y) + rho*r3/(8*sigma)*sqr(yDot));

        label n = 0;
        scalar rNew = 0.0;
        switch (SMDMethod_)
        {
		case method1:
                {
                    #include "NNTABSMDCalcMethod1.H"
                    break;
                }
                case method2:
                {
                    #include "NNTABSMDCalcMethod2.H"
                    break;
                }
        }

        if (rNew < r)
        {
                d = 2*rNew;
                y = 0;
                yDot = 0;
	        tc=0;
        }

    // update the nParticle count to conserve mass
    nParticle = semiMass/pow3(d);
    }

    if (nParticle > NP_)
    {
        addParcel = true;
        dChild = Factor1_*d;
        massChild = C_*semiMass;
        d = Factor2_*d;
    	// update the nParticle count to conserve mass
    	nParticle = ((1-C_)*semiMass)/pow3(d);
    }

    // Do not add child parcel
    return addParcel;
}
// ************************************************************************* //
