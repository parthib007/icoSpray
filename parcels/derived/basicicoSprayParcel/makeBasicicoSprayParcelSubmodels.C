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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "basicicoSprayCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeParcelForces.H" // thermo variant
//#include "makeThermoParcelTurbulenceForces.H" // add turbulence variant
#include "makeParcelDispersionModels.H"
#include "makeParcelTurbulenceDispersionModels.H" // add turbulence variant
#include "makeSprayParcelInjectionModels.H" // Spray variant
#include "makeParcelPatchInteractionModels.H"
//#include "makeParcelStochasticCollisionModels.H"
#include "makeParcelSurfaceFilmModels.H"
#include "makeSprayParcelStochasticCollisionModels.H" // Spray variant

// Spray
#include "DistortedSphereDragForce.H"
#include "makeSprayParcelAtomizationModels.H"
#include "makeSprayParcelBreakupModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicicoSprayCloud);

// Kinematic sub-models
makeParcelForces(basicicoSprayCloud);
//makeThermoParcelTurbulenceForces(basicicoSprayCloud);
makeParcelDispersionModels(basicicoSprayCloud);
makeParcelTurbulenceDispersionModels(basicicoSprayCloud);
makeSprayParcelInjectionModels(basicicoSprayCloud);
makeParcelPatchInteractionModels(basicicoSprayCloud);
//makeParcelStochasticCollisionModels(basicicoSprayCloud);
makeParcelSurfaceFilmModels(basicicoSprayCloud);
makeSprayParcelStochasticCollisionModels(basicicoSprayCloud);

// Spray sub-models
makeParticleForceModelType(DistortedSphereDragForce, basicicoSprayCloud);
makeSprayParcelAtomizationModels(basicicoSprayCloud);
makeSprayParcelBreakupModels(basicicoSprayCloud);


// ************************************************************************* //
