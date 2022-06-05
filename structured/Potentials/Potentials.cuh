#ifndef __POTENTIALS__
#define __POTENTIALS__

#include "CommonPotentials.cuh"
#include "CommonParameters.cuh"

#include "Bond1/Bond1.cuh"
#include "Bond2/Bond2.cuh"
#include "Bond3/Bond3.cuh"
#include "Bond4/Bond4.cuh"

//#include "UnBound/Steric.cuh"
//#include "UnBound/HarmonicConstraint.cuh"
#include "UnBound/LennardJones.cuh"
//#include "UnBound/WCA.cuh"
//#include "UnBound/LJ_WCA.cuh"
#include "UnBound/DebyeHuckel.cuh"
//#include "UnBound/Statistical.cuh"
#include "UnBound/KimHummer.cuh"
//#include "UnBound/Ravikumar.cuh"
#include "UnBound/Clashed.cuh"

#include "Surface/GenericSurface.cuh"
//#include "Surface/KaranicolasBrooksSurfacePotential.cuh"

#include "Bounds/SphericalShell.cuh"

#endif
