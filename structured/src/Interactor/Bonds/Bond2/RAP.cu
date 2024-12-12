#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond2/Bond2.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/RotationalAlignmentPotential.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct RAP_{

        struct ComputationalData{
            real4* pos;
            real4* dir;

            Box    box;
        };

        //Potential parameters
        struct StorageData{
        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.dir = pd->getDir(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;
            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real e;
            return e;
        }

        static inline __device__ ForceTorque forceTorque(int index_i, int index_j,
                                                         int currentParticleIndex,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){
            ForceTorque forceTorque;
            return forceTorque;
        }


    };

    using RAP = Bond2Torque_<RAP_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond2,RAP,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::RAP>
)
