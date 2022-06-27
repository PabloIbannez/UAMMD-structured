#ifndef MEASURES_BASIC_CUH
#define MEASURES_BASIC_CUH

#include "MeasuresTransforms.cuh"

namespace uammd{
namespace structured{
namespace Measures{
    
    real totalVirial(std::shared_ptr<ParticleGroup> pg,
                      cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 
        
        int N = pg->getNumberParticles();

        real* virial = pd->getVirial(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::totalVirial pp(virial);

        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        real tpp = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                  thrust::make_transform_iterator(pgIter, pp),
                                  thrust::make_transform_iterator(pgIter + N, pp),real(0.0));

        cudaStreamSynchronize(st);

        return tpp;
    }

    real totalVirial(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real v = totalVirial(pg,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    tensor3  totalStress(std::shared_ptr<ParticleGroup> pg,
                         cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 
        
        int N = pg->getNumberParticles();

        tensor3* stress = pd->getStress(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::totalStress pp(stress);

        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        tensor3 tpp = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                     thrust::make_transform_iterator(pgIter, pp),
                                     thrust::make_transform_iterator(pgIter + N, pp),(tensor3){0,0,0,
                                                                                               0,0,0,
                                                                                               0,0,0});
        
        cudaStreamSynchronize(st);

        return tpp;
    }
    
    tensor3 totalStress(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        tensor3 v = totalStress(pg,0);
        cudaDeviceSynchronize();
        return v;
    }

    tensor3  totalKineticPressure(std::shared_ptr<ParticleGroup> pg,
                                  cudaStream_t st){
        
        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 
        
        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::kineticPressure kp(mass,vel);

        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        tensor3 tkp = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                     thrust::make_transform_iterator(pgIter, kp),
                                     thrust::make_transform_iterator(pgIter + N, kp),(tensor3){0,0,0,
                                                                                               0,0,0,
                                                                                               0,0,0});

        cudaStreamSynchronize(st);
        
        return tkp;
    }

    tensor3 totalKineticPressure(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        tensor3 v = totalKineticPressure(pg,0);
        cudaDeviceSynchronize();
        return v;
    }

    real  totalKineticEnergy(std::shared_ptr<ParticleGroup> pg,
                             cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 
        
        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::kineticEnergy ke(mass,vel);

        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        real tke = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                  thrust::make_transform_iterator(pgIter, ke),
                                  thrust::make_transform_iterator(pgIter + N, ke),(real){0});
        
        cudaStreamSynchronize(st);

        return tke;
    }

    real totalKineticEnergy(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real v = totalKineticEnergy(pg,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real maxForce(std::shared_ptr<ParticleGroup> pg,
                  cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real4* force = pd->getForce(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::maxForce mF(force);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        real maxForce = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                       thrust::make_transform_iterator(pgIter, mF),
                                       thrust::make_transform_iterator(pgIter+N, mF),(real){-1.0},thrust::maximum<real>());
        
        cudaStreamSynchronize(st);

        return maxForce;
    }
    
    real maxForce(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real v = maxForce(pg,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real totalPotentialEnergy(std::shared_ptr<ParticleGroup> pg,
                              cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real* energy = pd->getEnergy(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::potentialEnergy pE(energy);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        real tEnergy = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                      thrust::make_transform_iterator(pgIter  , pE),
                                      thrust::make_transform_iterator(pgIter+N, pE),(real){0});
        
        cudaStreamSynchronize(st);

        return tEnergy;
    }
    
    real totalPotentialEnergy(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real v = totalPotentialEnergy(pg,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real3 totalForce(std::shared_ptr<ParticleGroup> pg,
                     cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real4* force = pd->getForce(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::totalForce tF(force);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        real3 tForce = make_real3(thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                       thrust::make_transform_iterator(pgIter, tF),
                       thrust::make_transform_iterator(pgIter+N, tF),(real4){0,0,0,0}));
        
        cudaStreamSynchronize(st);

        return tForce;
    }
    
    real3 totalForce(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real3 v = totalForce(pg,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real3 centroidPos(std::shared_ptr<ParticleGroup> pg,
                      cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::totalPos tP(pos);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 centroid = make_real3(thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                         thrust::make_transform_iterator(pgIter, tP),
                         thrust::make_transform_iterator(pgIter + N, tP),(real4){0,0,0,0})/N);
        
        cudaStreamSynchronize(st);

        return centroid;
    }
    
    real3 centroidPos(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real3 v = centroidPos(pg,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real totalCharge(std::shared_ptr<ParticleGroup> pg,
                     cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real* charge = pd->getCharge(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::totalCharge tQ(charge);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        real tCharge = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                     thrust::make_transform_iterator(pgIter, tQ),
                     thrust::make_transform_iterator(pgIter+N, tQ),(real){0});
        
        cudaStreamSynchronize(st);

        return tCharge;
    }
    
    real totalCharge(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real v = totalCharge(pg,0);
        cudaDeviceSynchronize();
        return v;
    }

    real totalMass(std::shared_ptr<ParticleGroup> pg,
                   cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::totalMass tM(mass);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);
        
        real tMass = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                     thrust::make_transform_iterator(pgIter, tM),
                     thrust::make_transform_iterator(pgIter+N, tM),(real){0});
        
        cudaStreamSynchronize(st);

        return tMass;
    }
    
    real totalMass(std::shared_ptr<ParticleGroup> pg){
        cudaDeviceSynchronize();
        real v = totalMass(pg,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real3 centerOfMassPos(std::shared_ptr<ParticleGroup> pg,
                          real   totalMass,
                          cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::weightedSum<real4> mWs(mass,pos);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 comp = make_real3(thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                thrust::make_transform_iterator(pgIter, mWs),
                                thrust::make_transform_iterator(pgIter + N, mWs),(real4){0,0,0,0})/totalMass);
        
        cudaStreamSynchronize(st);

        return comp;
    }
    
    real3 centerOfMassPos(std::shared_ptr<ParticleGroup> pg,
                          real totalMass){
        cudaDeviceSynchronize();
        real3 v = centerOfMassPos(pg,totalMass,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real3 centerOfMassVel(std::shared_ptr<ParticleGroup> pg,
                          real   totalMass,
                          cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::weightedSum<real3> mWs(mass,vel);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 comv = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                                    thrust::make_transform_iterator(pgIter, mWs),
                                    thrust::make_transform_iterator(pgIter + N, mWs),(real3){0,0,0})/totalMass;
        
        cudaStreamSynchronize(st);

        return comv;
    }
    
    real3 centerOfMassVel(std::shared_ptr<ParticleGroup> pg,
                          real totalMass){
        cudaDeviceSynchronize();
        real3 v = centerOfMassVel(pg,totalMass,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real3 angularMomentum(std::shared_ptr<ParticleGroup> pg,
                          real3 refp,
                          real3 refv,
                          cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::angularMomentum aM(mass,
                                               pos,
                                               vel,
                                               refp,refv);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 angM = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                     thrust::make_transform_iterator(pgIter, aM),
                     thrust::make_transform_iterator(pgIter + N, aM),(real3){0,0,0});
        
        cudaStreamSynchronize(st);

        return angM;
    }
    
    real3 angularMomentum(std::shared_ptr<ParticleGroup> pg,
                          real3 refp,
                          real3 refv){
        cudaDeviceSynchronize();
        real3 v = angularMomentum(pg,refp,refv,0);
        cudaDeviceSynchronize();
        return v;
    }

    tensor3 inertiaTensor(std::shared_ptr<ParticleGroup> pg,
                          real3 refp,
                          cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        int N = pg->getNumberParticles();
        
        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
        
        MeasuresTransforms::inertiaTensor I(mass,pos,refp);
        
        auto pgIter = pg->getIndexIterator(access::location::gpu);

        tensor3 inertia = thrust::reduce(thrust::cuda::par(sys->getTemporaryDeviceAllocator<char>()).on(st),
                          thrust::make_transform_iterator(pgIter, I),
                          thrust::make_transform_iterator(pgIter + N, I),tensor3(0));
        
        cudaStreamSynchronize(st);

        return inertia;
    }
    
    tensor3 inertiaTensor(std::shared_ptr<ParticleGroup> pg,
                          real3 refp,
                          real3 refv){
        cudaDeviceSynchronize();
        tensor3 v = inertiaTensor(pg,refp,0);
        cudaDeviceSynchronize();
        return v;
    }
    
    real3 angularVelocity(std::shared_ptr<ParticleGroup> pg,
                          real3  refp,
                          real3  refv,
                          cudaStream_t st){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem(); 

        tensor3 inertia = inertiaTensor(pg,refp,st);

        real m11 = inertia.xx;
        real m12 = inertia.xy;
        real m13 = inertia.xz;

        real m22 = inertia.yy;
        real m23 = inertia.yz;

        real m33 = inertia.zz;
        
        real D = m11*(m33*m22-m23*m23)-m12*(m33*m12-m23*m13)+m13*(m23*m12-m22*m13);

        real a11 = m33*m22-m23*m23;
        real a12 = m13*m23-m33*m12;  
        real a13 = m12*m23-m13*m22;  

        real a22 = m33*m11-m13*m13;
        real a23 = m12*m13-m11*m23;

        real a33 = m11*m22-m12*m12;

        D = (m11*a11)+(m12*a12)+(m13*a13); 

        tensor3 invInertia(a11,a12,a13,
                           a12,a22,a23,
                           a13,a23,a33);

        invInertia=invInertia/D;
        
        real3 angm = angularMomentum(pg,refp,refv,st);

        real3 angularVelocity;

        angularVelocity.x = invInertia.xx*angm.x+invInertia.xy*angm.y+invInertia.xz*angm.z;
        angularVelocity.y = invInertia.yx*angm.x+invInertia.yy*angm.y+invInertia.yz*angm.z;
        angularVelocity.z = invInertia.zx*angm.x+invInertia.zy*angm.y+invInertia.zz*angm.z;
        
        cudaStreamSynchronize(st);

        return angularVelocity;
    }
    
    real3 angularVelocity(std::shared_ptr<ParticleGroup> pg,
                          real3 refp,
                          real3 refv){
        cudaDeviceSynchronize();
        real3 v = angularVelocity(pg,refp,refv,0);
        cudaDeviceSynchronize();
        return v;
    }

}}}


#endif
