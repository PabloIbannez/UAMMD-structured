#include "UAMMDstructuredBase.cuh"
#include "utils/quaternion.cuh"
#include "Utils/Measures/MeasuresTransforms.cuh"

namespace uammd{
namespace structured{
namespace Measures{

    real meanDistance(std::shared_ptr<ParticleGroup> pg,
                      real3 point,
                      Box box,
                      cudaStream_t st = 0){

        auto pd  = pg->getParticleData();

        int N = pg->getNumberParticles();

        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::totalDistance r(pos,point,box);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real mr = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                 thrust::make_transform_iterator(pgIter, r),
                                 thrust::make_transform_iterator(pgIter + N, r),real(0.0));

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return mr/N;
    }

    real gyrationRadius(std::shared_ptr<ParticleGroup> pg,
                        real3 center,
                        Box box,
                        cudaStream_t st = 0){

        auto pd  = pg->getParticleData();

        int N = pg->getNumberParticles();

        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::totalSquareDistance r2(pos,center,box);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real mr = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                 thrust::make_transform_iterator(pgIter, r2),
                                 thrust::make_transform_iterator(pgIter + N, r2),real(0.0));

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return std::sqrt(mr/N);
    }

    real2 extremalDistances(std::shared_ptr<ParticleGroup> pg,
                            real3 point,
                            Box box,
                            cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::totalDistance r(pos,point,box);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real min = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                  thrust::make_transform_iterator(pgIter, r),
                                  thrust::make_transform_iterator(pgIter + N, r),(real){ INFINITY},thrust::minimum<real>());

        real max = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                  thrust::make_transform_iterator(pgIter, r),
                                  thrust::make_transform_iterator(pgIter + N, r),(real){-INFINITY},thrust::maximum<real>());

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return {min,max};
    }

    std::tuple<real2,real2,real2> extremalPositions(std::shared_ptr<ParticleGroup> pg,
                                                    Box box,
                                                    cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        MeasuresTransforms::posX px(pos,box);
        MeasuresTransforms::posY py(pos,box);
        MeasuresTransforms::posZ pz(pos,box);

        real minX = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                   thrust::make_transform_iterator(pgIter, px),
                                   thrust::make_transform_iterator(pgIter + N, px),(real){ INFINITY},thrust::minimum<real>());

        real maxX = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                   thrust::make_transform_iterator(pgIter, px),
                                   thrust::make_transform_iterator(pgIter + N, px),(real){-INFINITY},thrust::maximum<real>());

        real minY = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                   thrust::make_transform_iterator(pgIter, py),
                                   thrust::make_transform_iterator(pgIter + N, py),(real){ INFINITY},thrust::minimum<real>());

        real maxY = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                   thrust::make_transform_iterator(pgIter, py),
                                   thrust::make_transform_iterator(pgIter + N, py),(real){-INFINITY},thrust::maximum<real>());

        real minZ = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                   thrust::make_transform_iterator(pgIter, pz),
                                   thrust::make_transform_iterator(pgIter + N, pz),(real){ INFINITY},thrust::minimum<real>());

        real maxZ = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                   thrust::make_transform_iterator(pgIter, pz),
                                   thrust::make_transform_iterator(pgIter + N, pz),(real){-INFINITY},thrust::maximum<real>());

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return std::make_tuple<real2,real2,real2>({minX,maxX},
                                                  {minY,maxY},
                                                  {minZ,maxZ});
    }

    tensor3  totalKineticPressure(std::shared_ptr<ParticleGroup> pg,
                                  cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::kineticPressure kp(mass,vel);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        tensor3 tkp = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                     thrust::make_transform_iterator(pgIter, kp),
                                     thrust::make_transform_iterator(pgIter + N, kp),(tensor3){0,0,0,
                                                                                               0,0,0,
                                                                                               0,0,0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tkp;
    }

    real  totalKineticEnergy(std::shared_ptr<ParticleGroup> pg,
                                  cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::kineticEnergy ke(mass,vel);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real tke = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                  thrust::make_transform_iterator(pgIter, ke),
                                  thrust::make_transform_iterator(pgIter + N, ke),(real){0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tke;
    }

    real maxForce(std::shared_ptr<ParticleGroup> pg,
                  cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real4* force = pd->getForce(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::maxForce mF(force);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real maxForce = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                       thrust::make_transform_iterator(pgIter, mF),
                                       thrust::make_transform_iterator(pgIter+N, mF),(real){-INFINITY},thrust::maximum<real>());

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return maxForce;
    }

    real totalPotentialEnergy(std::shared_ptr<ParticleGroup> pg,
                              cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* energy = pd->getEnergy(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::potentialEnergy pE(energy);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real tEnergy = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                      thrust::make_transform_iterator(pgIter  , pE),
                                      thrust::make_transform_iterator(pgIter+N, pE),(real){0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tEnergy;
    }

    real totalLambdaDerivative(std::shared_ptr<ParticleGroup> pg,
                               cudaStream_t st = 0){

        auto pd  = pg->getParticleData();

        int N = pg->getNumberParticles();

        real* lambdaDerivative = pd->getLambdaDerivative(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::lambdaDerivative lD(lambdaDerivative);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real tLambdaDerivative = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                                thrust::make_transform_iterator(pgIter  , lD),
                                                thrust::make_transform_iterator(pgIter+N, lD),(real){0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tLambdaDerivative;

    }

    real3 totalForce(std::shared_ptr<ParticleGroup> pg,
                     cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real4* force = pd->getForce(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::totalForce tF(force);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 tForce = make_real3(thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                       thrust::make_transform_iterator(pgIter, tF),
                       thrust::make_transform_iterator(pgIter+N, tF),(real4){0,0,0,0}));

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tForce;
    }

    real3 centroidPos(std::shared_ptr<ParticleGroup> pg,
                      cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::totalPos tP(pos);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 centroid = make_real3(thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                    thrust::make_transform_iterator(pgIter, tP),
                                    thrust::make_transform_iterator(pgIter + N, tP),(real4){0,0,0,0})/N);

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return centroid;
    }

    real totalCharge(std::shared_ptr<ParticleGroup> pg,
                     cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* charge = pd->getCharge(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::totalCharge tQ(charge);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real tCharge = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                     thrust::make_transform_iterator(pgIter, tQ),
                     thrust::make_transform_iterator(pgIter+N, tQ),(real){0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tCharge;
    }

    real totalMass(std::shared_ptr<ParticleGroup> pg,
                   cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::totalMass tM(mass);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real tMass = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                     thrust::make_transform_iterator(pgIter, tM),
                     thrust::make_transform_iterator(pgIter+N, tM),(real){0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tMass;
    }

    real3 centerOfMassPos(std::shared_ptr<ParticleGroup> pg,
                          real totalMass,
                          cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::weightedSum<real4> mWs(mass,pos);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 comp = make_real3(thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                thrust::make_transform_iterator(pgIter, mWs),
                                thrust::make_transform_iterator(pgIter + N, mWs),(real4){0,0,0,0}));

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return comp/totalMass;
    }

    real3 centerOfMassVel(std::shared_ptr<ParticleGroup> pg,
                          real   totalMass,
                          cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::weightedSum<real3> mWs(mass,vel);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 comv = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                                    thrust::make_transform_iterator(pgIter, mWs),
                                    thrust::make_transform_iterator(pgIter + N, mWs),(real3){0,0,0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return comv/totalMass;
    }

    real3 angularMomentum(std::shared_ptr<ParticleGroup> pg,
                          real3 refp,
                          real3 refv,
                          Box box,
                          cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();
        real3* vel = pd->getVel(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::angularMomentum aM(mass,
                                               pos,
                                               vel,
                                               refp,refv,
                                               box);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 angM = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                     thrust::make_transform_iterator(pgIter, aM),
                     thrust::make_transform_iterator(pgIter + N, aM),(real3){0,0,0});

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return angM;
    }

    tensor3 inertiaTensor(std::shared_ptr<ParticleGroup> pg,
                          real3 refp,
                          Box box,
                          cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real* mass = pd->getMass(access::location::gpu, access::mode::read).raw();
        real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

        MeasuresTransforms::inertiaTensor I(mass,pos,refp,box);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        tensor3 inertia = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                          thrust::make_transform_iterator(pgIter, I),
                          thrust::make_transform_iterator(pgIter + N, I),tensor3(0));

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return inertia;
    }

    real3 angularVelocity(std::shared_ptr<ParticleGroup> pg,
                          real3  refp,
                          real3  refv,
                          Box    box,
                          cudaStream_t st = 0){

        auto pd  = pg->getParticleData();

        tensor3 inertia = inertiaTensor(pg,refp,box,st);

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

        real3 angm = angularMomentum(pg,refp,refv,box,st);

        real3 angularVelocity;

        angularVelocity.x = invInertia.xx*angm.x+invInertia.xy*angm.y+invInertia.xz*angm.z;
        angularVelocity.y = invInertia.yx*angm.x+invInertia.yy*angm.y+invInertia.yz*angm.z;
        angularVelocity.z = invInertia.zx*angm.x+invInertia.zy*angm.y+invInertia.zz*angm.z;

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return angularVelocity;
    }

      real3 totalMagnetization(std::shared_ptr<ParticleGroup> pg,
			      cudaStream_t st = 0){

        auto pd  = pg->getParticleData();


        int N = pg->getNumberParticles();

        real4* m_and_M = pd->getMagnetization(access::location::gpu, access::mode::read).raw();
	real4* dir = pd->getDir(access::location::gpu, access::mode::read).raw();
        MeasuresTransforms::magneticMoment_vec tM(dir, m_and_M);

        auto pgIter = pg->getIndexIterator(access::location::gpu);

        real3 tMagn = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
				     thrust::make_transform_iterator(pgIter, tM),
				     thrust::make_transform_iterator(pgIter+N, tM), real3());

        if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

        return tMagn;
      }

      real maxMagnetization(std::shared_ptr<ParticleGroup> pg,
              cudaStream_t st = 0){

          auto pd  = pg->getParticleData();


          int N = pg->getNumberParticles();

          real4* m_and_M = pd->getMagnetization(access::location::gpu, access::mode::read).raw();

          auto pgIter = pg->getIndexIterator(access::location::gpu);

          MeasuresTransforms::magneticMoment_mod tMmod(m_and_M);
          real tMagnMod = thrust::reduce(thrust::cuda::par(System::getTemporaryDeviceAllocator<char>()).on(st),
                  thrust::make_transform_iterator(pgIter, tMmod),
                  thrust::make_transform_iterator(pgIter+N, tMmod), 0.0);
          if(st==0){cudaDeviceSynchronize();} else {cudaStreamSynchronize(st);}

          return tMagnMod;
      }
}}}
