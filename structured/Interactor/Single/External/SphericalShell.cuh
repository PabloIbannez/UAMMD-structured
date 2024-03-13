#ifndef __EXTERNAL_SPHERICAL_SHELL_POT__
#define __EXTERNAL_SPHERICAL_SHELL_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct SphericalShell_{

        //Computational data
        struct ComputationalData{

            real4* pos;
            real*  radius;

            real3 shellCenter;
            real  shellRadius;

            real shellEpsilon;
            real shellSigma;
        };

        struct StorageData{
            real3 shellCenter;
            real  shellRadius;

            real shellEpsilon;
            real shellSigma;

            real minShellRadius;
            real maxShellRadius;

            real radiusVelocity;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.radius = pd->getRadius(access::location::gpu, access::mode::read).raw();

            computational.shellCenter = storage.shellCenter;

            computational.shellEpsilon = storage.shellEpsilon;
            computational.shellSigma   = storage.shellSigma;

            //Computing the shell radius

            double simulationTime = gd->getFundamental()->getSimulationTime();

            computational.shellRadius = storage.shellRadius + storage.radiusVelocity*simulationTime;

            if(computational.shellRadius < storage.minShellRadius){
                computational.shellRadius = storage.minShellRadius;
            }

            if(computational.shellRadius > storage.maxShellRadius){
                computational.shellRadius = storage.maxShellRadius;
            }

            System::log<System::DEBUG>("[SphericalShell] Current shell radius: %f", computational.shellRadius);

            return computational;
        }

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.shellCenter = data.getParameter<real3>("shellCenter");
            storage.shellRadius = data.getParameter<real>("shellRadius");

            storage.shellEpsilon = data.getParameter<real>("shellEpsilon",1.0);
            storage.shellSigma   = data.getParameter<real>("shellSigma",1.0);

            storage.minShellRadius = data.getParameter<real>("minShellRadius",0.0);
            storage.maxShellRadius = data.getParameter<real>("maxShellRadius",INFINITY);

            storage.radiusVelocity = data.getParameter<real>("radiusVelocity",0.0);

            if(storage.minShellRadius > storage.maxShellRadius){
                System::log<System::CRITICAL>("[SphericalShell] Minimum shell radius (%f) is greater than maximum shell radius (%f).",
                                              storage.minShellRadius, storage.maxShellRadius);
            }

            if(storage.shellRadius < storage.minShellRadius){
                System::log<System::CRITICAL>("[SphericalShell] Shell radius (%f) is smaller than minimum shell radius (%f).",
                                              storage.shellRadius, storage.minShellRadius);
            }

            if(storage.shellRadius > storage.maxShellRadius){
                System::log<System::CRITICAL>("[SphericalShell] Shell radius (%f) is greater than maximum shell radius (%f).",
                                              storage.shellRadius, storage.maxShellRadius);
            }

            System::log<System::MESSAGE>("[SphericalShell] Shell center: (%f, %f, %f)",
                                          storage.shellCenter.x, storage.shellCenter.y, storage.shellCenter.z);
            System::log<System::MESSAGE>("[SphericalShell] Shell radius: %f", storage.shellRadius);
            System::log<System::MESSAGE>("[SphericalShell] Shell epsilon: %f", storage.shellEpsilon);
            System::log<System::MESSAGE>("[SphericalShell] Shell sigma: %f", storage.shellSigma);

            System::log<System::MESSAGE>("[SphericalShell] Minimum shell radius: %f", storage.minShellRadius);
            System::log<System::MESSAGE>("[SphericalShell] Maximum shell radius: %f", storage.maxShellRadius);

            System::log<System::MESSAGE>("[SphericalShell] Radius velocity: %f", storage.radiusVelocity);

            return storage;
        }

        static inline __device__ real energy(int index_i,const ComputationalData& computational){
            const real3 rij = make_real3(computational.pos[index_i]) - computational.shellCenter;

            real r2 = sqrt(dot(rij,rij)) - computational.shellRadius;
                 r2 = r2*r2;

            //rij,r2,epsilon,sigma
            const real sigma = computational.shellSigma+computational.radius[index_i];
            real e = BasicPotentials::Steric::Steric12::energy(rij,r2,computational.shellEpsilon,sigma);

            return e;
        }


        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
            const real3 rij = make_real3(computational.pos[index_i]) - computational.shellCenter;

            const real r = sqrt(dot(rij,rij));
                  real r2 = r - computational.shellRadius;
                       r2 = r2*r2;

            const real sinvr2  = computational.shellSigma*computational.shellSigma/r2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;
            const real sinvr12 = sinvr6*sinvr6;

            real fmod = real(12.0)*computational.shellEpsilon*sinvr12/(r-computational.shellRadius);
            real3 f;
            if(r <= std::numeric_limits<real>::epsilon()){
                f = make_real3(0.0);
            } else {
                f = fmod*rij/r;
            }

            return f;
        }

    };

    using SphericalShell = External_<SphericalShell_>;

}}}}

#endif
