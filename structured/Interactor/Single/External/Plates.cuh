#ifndef __EXTERNAL_PLATES_POT__
#define __EXTERNAL_PLATES_POT__

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

    struct Plates_{

        //Computational data
        struct ComputationalData{

            real4* pos;
            real*  radius;

            real topPlatePos;
            real bottomPlatePos;

            real platesEpsilon;
            real platesSigma;
        };

        struct StorageData{

            real compressionVelocity;

            real platesSeparation;

            real minPlatesSeparation;
            real maxPlatesSeparation;

            real platesEpsilon;
            real platesSigma;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.radius = pd->getRadius(access::location::gpu, access::mode::read).raw();

            computational.platesEpsilon = storage.platesEpsilon;
            computational.platesSigma   = storage.platesSigma;

            //Computing the top and bottom plate positions

            double simulationTime = gd->getFundamental()->getSimulationTime();

            computational.topPlatePos    = double( storage.platesSeparation/2.0) + double(storage.compressionVelocity)*simulationTime;
            computational.bottomPlatePos = double(-storage.platesSeparation/2.0) - double(storage.compressionVelocity)*simulationTime;

            real currentPlatesSeparation = computational.topPlatePos - computational.bottomPlatePos;

            if(currentPlatesSeparation < storage.minPlatesSeparation){
                computational.topPlatePos    =  storage.minPlatesSeparation/2.0;
                computational.bottomPlatePos = -storage.minPlatesSeparation/2.0;
            }

            if(currentPlatesSeparation > storage.maxPlatesSeparation){
                computational.topPlatePos    =  storage.maxPlatesSeparation/2.0;
                computational.bottomPlatePos = -storage.maxPlatesSeparation/2.0;
            }

            System::log<System::DEBUG>("[Plates] Current plates separation: %f", currentPlatesSeparation);

            return computational;
        }

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.platesSeparation = data.getParameter<real>("platesSeparation");

            storage.platesEpsilon = data.getParameter<real>("platesEpsilon", 1.0);
            storage.platesSigma   = data.getParameter<real>("platesSigma",   1.0);

            storage.compressionVelocity = data.getParameter<real>("compressionVelocity", 0.0);

            storage.minPlatesSeparation = data.getParameter<real>("minPlatesSeparation", 0.0);
            storage.maxPlatesSeparation = data.getParameter<real>("maxPlatesSeparation", INFINITY);

            if(storage.minPlatesSeparation > storage.maxPlatesSeparation){
                System::log<System::CRITICAL>("[Plates] Minimum plates separation is greater than maximum plates separation.");
            }

            if(storage.platesSeparation < storage.minPlatesSeparation){
                System::log<System::CRITICAL>("[Plates] Initial plates separation is smaller than minimum plates separation.");
            }

            if(storage.platesSeparation > storage.maxPlatesSeparation){
                System::log<System::CRITICAL>("[Plates] Initial plates separation is greater than maximum plates separation.");
            }

            System::log<System::MESSAGE>("[Plates] Initial plates separation: %f", storage.platesSeparation);

            System::log<System::MESSAGE>("[Plates] Minimum plates separation: %f", storage.minPlatesSeparation);
            System::log<System::MESSAGE>("[Plates] Maximum plates separation: %f", storage.maxPlatesSeparation);

            System::log<System::MESSAGE>("[Plates] Compression velocity: %f", storage.compressionVelocity);

            System::log<System::MESSAGE>("[Plates] Plates epsilon: %f", storage.platesEpsilon);
            System::log<System::MESSAGE>("[Plates] Plates sigma:   %f", storage.platesSigma);

            return storage;
        }

        static inline __device__ real energy(int index_i,const ComputationalData& computational){
            const real dzTop    = computational.topPlatePos    - computational.pos[index_i].z;
            const real dzBottom = computational.bottomPlatePos - computational.pos[index_i].z;

            const real sigma2invdzTop2    = computational.platesSigma*computational.platesSigma/(dzTop*dzTop);
            const real sigma2invdzBottom2 = computational.platesSigma*computational.platesSigma/(dzBottom*dzBottom);

            const real sigma2invdzTop6  = sigma2invdzTop2*sigma2invdzTop2*sigma2invdzTop2;
            const real sigma2invdzBot6  = sigma2invdzBottom2*sigma2invdzBottom2*sigma2invdzBottom2;

            const real sigma2invdzTop12 = sigma2invdzTop6*sigma2invdzTop6;
            const real sigma2invdzBot12 = sigma2invdzBot6*sigma2invdzBot6;

            //E = eps*(sigma/dzTop)^12 + eps*(sigma/dzBottom)^12
            real e = computational.platesEpsilon*(sigma2invdzTop12 + sigma2invdzBot12);

            return e;
        }


        static inline __device__ real3 force(int index_i,const ComputationalData& computational){
            const real dzTop    = computational.topPlatePos    - computational.pos[index_i].z;
            const real dzBottom = computational.bottomPlatePos - computational.pos[index_i].z;

            const real sigma2invdzTop2    = computational.platesSigma*computational.platesSigma/(dzTop*dzTop);
            const real sigma2invdzBottom2 = computational.platesSigma*computational.platesSigma/(dzBottom*dzBottom);

            const real sigma2invdzTop6  = sigma2invdzTop2*sigma2invdzTop2*sigma2invdzTop2;
            const real sigma2invdzBot6  = sigma2invdzBottom2*sigma2invdzBottom2*sigma2invdzBottom2;

            const real sigma2invdzTop12 = sigma2invdzTop6*sigma2invdzTop6;
            const real sigma2invdzBot12 = sigma2invdzBot6*sigma2invdzBot6;

            const real3 fTop    = {real(0.0),real(0.0),-computational.platesEpsilon*sigma2invdzTop12/dzTop};
            const real3 fBottom = {real(0.0),real(0.0),-computational.platesEpsilon*sigma2invdzBot12/dzBottom};

            return fTop + fBottom;
        }

    };

    using Plates         = External_<Plates_>;

}}}}

#endif
