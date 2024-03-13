#ifndef __INTEGRATOR_UTILS__
#define __INTEGRATOR_UTILS__

#include "Integrator/Integrator.cuh"

namespace uammd{
namespace structured{
namespace IntegratorUtils{

    void generateVelocity(shared_ptr<ParticleGroup> pg,
                          real kBT,
                          uint seed,
                          cudaStream_t stream){

            auto pd = pg->getParticleData();

            int N = pg->getNumberParticles();

            if(pd->isVelAllocated()){
                System::log<System::WARNING>("[GenerateVelocity] Velocity overwritten!");
            }

            auto groupIterator = pg->getIndexIterator(access::location::gpu);

            auto mass = pd->getMass(access::location::gpu, access::mode::read);

            auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
            auto vel = pd->getVel(access::location::gpu, access::mode::readwrite);

            uint seed_v = seed;

            real* mass_ptr    = mass.raw();

            real4* pos_ptr     = pos.raw();
            real3* vel_ptr     = vel.raw();

            real sigma_kBT = sqrt(kBT);
            thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                    [=] __host__ __device__ (int index){const real  sigma = sigma_kBT/sqrt(mass_ptr[index]);
                    Saru rng_s(index, seed_v);
                    const real3 noise = make_real3(rng_s.gf(0, sigma), rng_s.gf(0, sigma).x);
                    vel_ptr[index]=noise;});

            CudaSafeCall(cudaStreamSynchronize(stream));

        }

        void loadFrictionConstant(shared_ptr<ParticleGroup> pg,
                                  real frictionConstant){

            auto pd = pg->getParticleData();

            System::log<System::MESSAGE>("[IntegratorUtils] Loading friction constant (%f) to all particles.", frictionConstant);

            auto fricConst = pd->getFrictionConstant(access::location::cpu,
                                                     access::mode::write);

            auto id = pd->getId(access::location::cpu,
                                access::mode::read);

            auto groupIndex  = pg->getIndexIterator(access::location::cpu);
            auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

            fori(0,pg->getNumberParticles()){
                int index = sortedIndex[id[groupIndex[i]]];
                fricConst[i] = frictionConstant;
            }
        }

        void loadMobility(shared_ptr<ParticleGroup> pg,
                          real viscosity){

            auto pd = pg->getParticleData();

            System::log<System::MESSAGE>("[IntegratorUtils] Loading viscosity (%f) to all particles. Computing translational and rotational self mobility.", viscosity);

            auto radius = pd->getRadius(access::location::cpu, access::mode::read);

            auto mt = pd->getTranslationalSelfDiffusion(access::location::cpu, access::mode::write);
            auto mr = pd->getRotationalSelfDiffusion(access::location::cpu, access::mode::write);

            auto id = pd->getId(access::location::cpu,
                               access::mode::read);

            auto groupIndex  = pg->getIndexIterator(access::location::cpu);
            auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

            fori(0,pg->getNumberParticles()){
                int index = sortedIndex[id[groupIndex[i]]];
                mt[index] = real(1.0)/(real(6.0)*M_PI*viscosity*radius[index]);
                mr[index] = real(1.0)/(real(8.0)*M_PI*viscosity*radius[index]*radius[index]*radius[index]);
            }

        }

}}}

#endif
