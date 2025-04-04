#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Patches/NonBondedPatches/NonBondedPatches.cuh"
#include "Interactor/Patches/PatchesFactory.cuh"

#include "Interactor/BasicPotentials/DistanceSwitch.cuh"
#include "Interactor/BasicParameters/Pair/DistanceSwitch.cuh"
#include "Utils/ParameterHandler/PairParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBondedPatches{

    struct DistanceSwitchCosine_{

        using ParametersType        = typename BasicParameters::Pairs::DistanceSwitchCosine;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        ///////////////////////////

        //Computational data
        struct ComputationalData{

            real4* patchesPos;
            real4* patchesVector;

            Box box;

            ParametersPairsIterator paramPairIterator;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> swtParam;

            real cutOff;
        };

        static StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          std::shared_ptr<GlobalData>    patchesGd,
                                          std::shared_ptr<ParticleGroup> patchesPg,
                                          DataEntry& data){

            StorageData storage;

            storage.swtParam = std::make_shared<ParameterPairsHandler>(patchesGd,patchesPg,
                                                                       data);

            /////////////////////////////////////////////////////////

            auto pairsParam = storage.swtParam->getPairParameters();

            real cutOff = 0.0;
            for(auto p : pairsParam){
                cutOff = std::max(cutOff, p.second.rc);
            }

            System::log<System::MESSAGE>("[DistanceSwitchCosine] cutOff: %f" ,cutOff);

            storage.cutOff = cutOff;

            return storage;
        }

        static ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                      std::shared_ptr<ParticleGroup> pg,
                                                      std::shared_ptr<GlobalData>    patchesGd,
                                                      std::shared_ptr<ParticleGroup> patchesPg,
                                                      const StorageData&  storage,
                                                      const Computables& comp,
                                                      const cudaStream_t& st){

            ComputationalData computational;

            computational.box               = patchesGd->getEnsemble()->getBox();
            computational.paramPairIterator = storage.swtParam->getPairIterator();

            std::shared_ptr<ParticleData> patchesPd = patchesPg->getParticleData();

            computational.patchesPos    = patchesPd->getPos(access::location::gpu, access::mode::read).raw();
            computational.patchesVector = patchesPd->getPatchVector(access::location::gpu, access::mode::read).raw();

            return computational;
        }

        //index_i is the current particle
        static inline __device__ EnergyForceTorque energyForceTorque(const int& index_i,const int& index_j,
                                                                     const ComputationalData& computational){

            const real4 patch_i_pos = computational.patchesPos[index_i];
            const real4 patch_j_pos = computational.patchesPos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(patch_j_pos)-make_real3(patch_i_pos));

            const auto param = computational.paramPairIterator(index_i,index_j);

            const real E       = param.E;
            const real rc      = param.rc;
            const real K       = param.K;

            const real r2 = dot(rij, rij);

            real4 fe = make_real4(0.0);
            real3 t  = make_real3(0.0);

            if(r2<=rc*rc){
                fe = -E*BasicPotentials::DistanceSwitchCosine::forceEnergy(rij,r2,rc,K);
                t  =  cross(make_real3(computational.patchesVector[index_i]),make_real3(fe));
            }

            EnergyForceTorque eFrcTrq;

            eFrcTrq.energy = fe.w;
            eFrcTrq.force  = make_real4(make_real3(fe),0.0);
            eFrcTrq.torque = make_real4(t,0.0);


            return eFrcTrq;

        }

    };

    using DistanceSwitchCosine = NonBondedDynamicallyBondedPatchyParticles_<DistanceSwitchCosine_>;

}}}}

REGISTER_NONBONDED_PATCHES_INTERACTOR(
    NonBondedPatches,DistanceSwitchCosine,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBondedPatches::DistanceSwitchCosine>
)
