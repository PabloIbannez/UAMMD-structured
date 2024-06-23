#include "EulerMaruyamaRigidBody.cuh"

REGISTER_INTEGRATOR(
    Brownian,EulerMaruyamaRigidBody,
    uammd::structured::Integrator::NVT::Brownian::EulerMaruyamaRigidBody
)

//REGISTER_INTEGRATOR(
//    Brownian,EulerMaruyamaRigidBodyPatchesState,
//    uammd::structured::Integrator::NVT::Brownian::EulerMaruyamaRigidBodyPatchesState
//)

