#ifndef __NVT_INCLUDERS_CUH__
#define __NVT_INCLUDERS_CUH__

#include "BDHIDoublyPeriodic/DPStokes.cuh"
#include "BDHIDoublyPeriodic/Quasi2D.cuh"
#include "BDHIOpenBoundary/Cholesky.cuh"
#include "BDHIOpenBoundary/Lanczos.cuh"
#include "BDHITriplyPeriodic/FluctuatingImmersedBoundary.cuh"
#include "BDHITriplyPeriodic/ForceCouplingMethod.cuh"
#include "BDHITriplyPeriodic/PositivelySplitEwald.cuh"
#include "Brownian/EulerMaruyama.cuh"
#include "Brownian/EulerMaruyamaRigidBody.cuh"
#include "DPD/Verlet.cuh"
#include "FluctuatingHydrodynamics/CompressibleInertialCoupling.cuh"
#include "FluctuatingHydrodynamics/IncompressibleInertialCoupling.cuh"
#include "Langevin/BBK.cuh"
#include "Langevin/GJF.cuh"
#include "Magnetic/Magnetic.cuh"
#include "Magnetic/Brownian.cuh"
#include "Magnetic/Fixed.cuh"
#include "Magnetic/ForceCouplingMethod.cuh"
#include "SPH/Verlet.cuh"

#endif //__NVT_INCLUDERS_CUH__
