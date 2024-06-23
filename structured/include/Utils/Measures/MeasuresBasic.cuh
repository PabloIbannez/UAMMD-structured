#pragma once

#include "utils/quaternion.cuh"
#include "MeasuresTransforms.cuh"
#include <memory>
#include <tuple>

namespace uammd {
namespace structured {
namespace Measures {

real totalCharge(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);
real totalMass(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);

real3 totalMagnetization(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);
real maxMagnetization(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);

real3 centroidPos(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);
real3 centerOfMassPos(std::shared_ptr<ParticleGroup> pg, real totalMass, cudaStream_t st = 0);
real3 centerOfMassVel(std::shared_ptr<ParticleGroup> pg, real totalMass, cudaStream_t st = 0);

real  gyrationRadius(std::shared_ptr<ParticleGroup> pg, real3 center, Box box, cudaStream_t st = 0);
tensor3 inertiaTensor(std::shared_ptr<ParticleGroup> pg, real3 refp, Box box, cudaStream_t st = 0);

real  meanDistance  (std::shared_ptr<ParticleGroup> pg, real3 point,  Box box, cudaStream_t st = 0);
real2 extremalDistances(std::shared_ptr<ParticleGroup> pg, real3 point, Box box, cudaStream_t st = 0);
std::tuple<real2,real2,real2> extremalPositions(std::shared_ptr<ParticleGroup> pg, Box box, cudaStream_t st = 0);

real maxForce(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);
real3 totalForce(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);

real totalKineticEnergy(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);
real totalPotentialEnergy(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);

tensor3 totalKineticPressure(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);

real totalLambdaDerivative(std::shared_ptr<ParticleGroup> pg, cudaStream_t st = 0);

real3 angularVelocity(std::shared_ptr<ParticleGroup> pg, real3 refp, real3 refv, Box box, cudaStream_t st = 0);
real3 angularMomentum(std::shared_ptr<ParticleGroup> pg, real3 refp, real3 refv, Box box, cudaStream_t st = 0);
}}}
