#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "NonBondedPatches.cuh"
#include "Interactor/Patches/PatchesFactory.cuh"

#include "Interactor/BasicPotentials.cuh"
#include "Interactor/BasicParameters.cuh"
#include "Utils/ParameterHandler/PairParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBondedPatches{

    template<typename potential>
    struct Helix_{

        using ParametersType        = typename BasicParameters::Pairs::Helix<potential>;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        ///////////////////////////

        //Computational data
        struct ComputationalData{

            real4* dir;
            int*   patchesParentIndex;

            real4* patchesPos;
            real4* patchesVector;

            Box box;

            int startType;
            int endType;

            ParametersPairsIterator paramPairIterator;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> helixParam;

            int startType;
            int endType;

            real cutOff;

        };

        template <typename ParameterType>
        static void checkTypesAndPatches(std::shared_ptr<GlobalData>    gd,
                                         std::shared_ptr<ParticleGroup> pg,
                                         std::shared_ptr<GlobalData>    patchesGd,
                                         std::shared_ptr<ParticleGroup> patchesPg,
                                         std::string startTypeName,
                                         std::string endTypeName,
                                         const ParameterType& storage){
            {
                std::shared_ptr<ParticleData> patchesPd = patchesPg->getParticleData();

                auto patchesPos     = patchesPd->getPos(access::location::cpu, access::mode::read);
                auto patchesPosIter = patchesPg->getPropertyIterator(patchesPos,access::location::cpu);

                //Check all types are startType or endType
                for(int i=0;i<patchesPg->getNumberParticles();i++){
                    int tpy = int(patchesPosIter[i].w);
                    if((tpy != storage.startType) and (tpy != storage.endType)){
                        System::log<System::CRITICAL>("[Helix] The type of the patches has to be startType (%s) or endType (%s).",
                                                      startTypeName.c_str(),endTypeName.c_str());
                    }
                }
            }

            //Check each particles has got only two patches,
            //one of type start and the other one of type end
            {
                std::map<int,bool> hasStype;
                std::map<int,bool> hasEtype;

                std::shared_ptr<ParticleData> patchesPd = patchesPg->getParticleData();

                auto patchesPos          = patchesPd->getPos(access::location::cpu, access::mode::read);
                auto patchesParentId     = patchesPd->getModelId(access::location::cpu, access::mode::read);

                auto patchesPosIter      = patchesPg->getPropertyIterator(patchesPos,access::location::cpu);
                auto patchesParentIdIter = patchesPg->getPropertyIterator(patchesParentId,access::location::cpu);

                for(int i=0;i<patchesPg->getNumberParticles();i++){
                    int pId = patchesParentIdIter[i];

                    hasStype[pId]=false;
                    hasEtype[pId]=false;
                }

                for(int i=0;i<patchesPg->getNumberParticles();i++){

                    int tpy = int(patchesPosIter[i].w);
                    int pId = patchesParentIdIter[i];

                    if( tpy == storage.startType){
                        if(hasStype[pId]==false){
                            hasStype[pId]=true;
                        } else {
                            System::log<System::CRITICAL>("[Helix] S type has been added before for particle %i.",pId);
                        }
                    } else {
                        if(hasEtype[pId]==false){
                            hasEtype[pId]=true;
                        } else {
                            System::log<System::CRITICAL>("[Helix] E type has been added before for particle %i.",pId);
                        }
                    }
                }

                for(int i=0;i<patchesPg->getNumberParticles();i++){
                    int pId = patchesParentIdIter[i];
                    if(!hasStype[pId]){
                        System::log<System::CRITICAL>("[Helix] S type has not been added for particle %i.",pId);
                    }
                    if(!hasEtype[pId]){
                        System::log<System::CRITICAL>("[Helix] E type has not been added for particle %i.",pId);
                    }
                }
            }
        }

        static StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          std::shared_ptr<GlobalData>    patchesGd,
                                          std::shared_ptr<ParticleGroup> patchesPg,
                                          DataEntry& data){

            StorageData storage;

            std::string startTypeName = data.getParameter<std::string>("startType","S");
            std::string endTypeName   = data.getParameter<std::string>("endType","E");

            System::log<System::MESSAGE>("[Helix] Added parameter startType: %s",startTypeName.c_str());
            System::log<System::MESSAGE>("[Helix] Added parameter endType: %s",endTypeName.c_str());

            storage.helixParam = std::make_shared<ParameterPairsHandler>(patchesGd,patchesPg,
                                                                         data);

            ////////////////////////////////////////////////////////////////////////////

            auto pairsParam = storage.helixParam->getPairParameters();

            real cutOff = 0.0;
            for(auto p : pairsParam){
                cutOff = std::max(cutOff, p.second.helixParams.dstParams.rc);
            }

            System::log<System::MESSAGE>("[Helix] cutOff: %f" ,cutOff);

            storage.cutOff = cutOff;

            ////////////////////////////////////////////////////////////////////////////

            //Check types

            auto types = patchesGd->getTypes();

            storage.startType = types->getTypeId(startTypeName);
            storage.endType   = types->getTypeId(endTypeName);

            ///////////////////////////////////////////////////////////

            checkTypesAndPatches<StorageData>(gd,pg,
                                              patchesGd,patchesPg,
                                              startTypeName,endTypeName,
                                              storage);

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

            computational.box = patchesGd->getEnsemble()->getBox();
            computational.paramPairIterator = storage.helixParam->getPairIterator();

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.dir                = pd->getDir(access::location::gpu, access::mode::read).raw();

            std::shared_ptr<ParticleData> patchesPd = patchesPg->getParticleData();
            computational.patchesParentIndex  = patchesPd->getParentIndex(access::location::gpu, access::mode::read).raw();

            computational.patchesPos          = patchesPd->getPos(access::location::gpu, access::mode::read).raw();
            computational.patchesVector       = patchesPd->getPatchVector(access::location::gpu, access::mode::read).raw();

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////

            computational.startType = storage.startType;
            computational.endType   = storage.endType;

            return computational;
        }

        static inline __device__ int setUpBond(const int& index_i,const int& index_j,
                                               const real4* __restrict__ patchesPos,
                                               const real4* __restrict__ patchesVector,
                                               const int*   __restrict__ patchesParentIndex,
                                               const real4* __restrict__ dir,
                                               const Box& box,
                                               const int startType,const int endType,
                                               Quat& q_s,Quat& q_e,
                                               real3& lps, real3& drse){
            // This function set the values of qs,qe,lps and drse
            // and return the index of the start patch

            // index_i (which is the index of the particle we are working on) is always considered as the be the start patch
            // we only have to determine the trully start patch. If index_i is the trully start patch then nothing is done
            // otherwise we have to transpose the R_H matrix.

            int index_s = -1;

            const real4 patch_i_pos = patchesPos[index_i];
            const real4 patch_j_pos = patchesPos[index_j];

            const int tpy_i = int(patch_i_pos.w);
            const int tpy_j = int(patch_j_pos.w);

            if ( tpy_i == startType and tpy_j == endType ){
                index_s = index_i;
            } else {
                index_s = index_j;
            }

            q_s = dir[patchesParentIndex[index_i]];
            q_e = dir[patchesParentIndex[index_j]];

            lps  = make_real3(patchesVector[index_i]);

            drse = box.apply_pbc(make_real3(patch_j_pos)-make_real3(patch_i_pos));

            //printf("startType %i, endType %i, \\
            //        index_i: %i, index_j: %i, \\
            //        tpy_i: %i, tpy_j: %i, index_s: %i, lps: %f %f %f, drse: %f %f %f\n",
            //        startType,endType,
            //        index_i,index_j,tpy_i,tpy_j,index_s,lps.x,lps.y,lps.z,drse.x,drse.y,drse.z);

            return index_s;
        }

        //index_i is the current particle
        static inline __device__ EnergyForceTorque energyForceTorque(const int& index_i,const int& index_j,
                                                                     const ComputationalData& computational){


            Quat qs,qe;
            real3 lps,drse;

            int index_s = setUpBond(index_i,index_j,
                                    computational.patchesPos,
                                    computational.patchesVector,
                                    computational.patchesParentIndex,
                                    computational.dir,
                                    computational.box,
                                    computational.startType,computational.endType,
                                    qs,qe,lps,drse);

            const auto param = computational.paramPairIterator(index_i,index_j);

            real3 es_x = qs.getVx();

            tensor3 R_H = param.R_H;
            if (index_s == index_j){
                R_H = R_H.transpose();
            }

            //At this point we have the following information:
            // - q_s, q_e: the orientation of the patches
            // - lps: the vector from the center of the particle to the patch for the start particle
            // - drse: the vector from the patch of the start particle to the patch of the end particle
            // - R_H: the rotation matrix from the local frame of the start particle to the local frame of the end particle

            ///////////////////////////////////////////////////

            const real Eb = param.Eb;

            const typename potential::params p = param.helixParams;

            ///////////////////////////////////////////////////

            EnergyForceTorque eft;

            eft.energy = BasicPotentials::Helix::Dynamic::energy<potential>(drse,qs,qe,R_H,
                                                                            p.dstParams,
                                                                            p.thetaParams,
                                                                            p.phiParams,
                                                                            Eb);

            eft.force  = make_real4(BasicPotentials::Helix::Dynamic::force<potential>
                                   (drse,qs,qe,R_H,
                                    p.dstParams,
                                    p.thetaParams,
                                    p.phiParams,
                                    Eb),0.0);

            eft.torque = make_real4(BasicPotentials::Helix::Dynamic::torque<potential>
                                   (drse,lps,qs,qe,R_H,
                                    p.dstParams,
                                    p.thetaParams,
                                    p.phiParams,
                                    Eb),0.0);

            return eft;
        }
    };

    template<typename potential>
    struct Helix2States_{

        using ParametersType        = typename BasicParameters::Pairs::Helix2States<potential>; // !!!
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        ///////////////////////////

        //Computational data
        struct ComputationalData{

            real4* dir;
            int*   patchesParentIndex;

            real4* patchesPos;
            real4* patchesVector;

            int4* patchesState;

            Box box;

            int startType;
            int endType;

            ParametersPairsIterator paramPairIterator;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> helixParam;

            int startType;
            int endType;

            real cutOff;

        };

        static StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          std::shared_ptr<GlobalData>    patchesGd,
                                          std::shared_ptr<ParticleGroup> patchesPg,
                                          DataEntry& data){

            StorageData storage;

            std::string startTypeName = data.getParameter<std::string>("startType","S");
            std::string endTypeName   = data.getParameter<std::string>("endType","E");

            System::log<System::MESSAGE>("[Helix] Added parameter startType: %s",startTypeName.c_str());
            System::log<System::MESSAGE>("[Helix] Added parameter endType: %s",endTypeName.c_str());

            storage.helixParam = std::make_shared<ParameterPairsHandler>(patchesGd,patchesPg,
                                                                         data);

            ////////////////////////////////////////////////////////////////////////////

            auto pairsParam = storage.helixParam->getPairParameters();

            real cutOff = 0.0;
            for(auto p : pairsParam){
                cutOff = std::max(cutOff, p.second.helixParams0.dstParams.rc);
                cutOff = std::max(cutOff, p.second.helixParams1.dstParams.rc);
            }

            System::log<System::MESSAGE>("[Helix] cutOff: %f" ,cutOff);

            storage.cutOff = cutOff;

            ////////////////////////////////////////////////////////////////////////////

            //Check types

            auto types = patchesGd->getTypes();

            storage.startType = types->getTypeId(startTypeName);
            storage.endType   = types->getTypeId(endTypeName);

            ///////////////////////////////////////////////////////////

            Helix_<potential>::template checkTypesAndPatches<StorageData>(gd,pg,
                                                                          patchesGd,patchesPg,
                                                                          startTypeName,endTypeName,
                                                                          storage);

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

            computational.box = patchesGd->getEnsemble()->getBox();
            computational.paramPairIterator = storage.helixParam->getPairIterator();

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.dir                = pd->getDir(access::location::gpu, access::mode::read).raw();

            std::shared_ptr<ParticleData> patchesPd = patchesPg->getParticleData();
            computational.patchesParentIndex  = patchesPd->getParentIndex(access::location::gpu, access::mode::read).raw();

            computational.patchesPos          = patchesPd->getPos(access::location::gpu, access::mode::read).raw();
            computational.patchesVector       = patchesPd->getPatchVector(access::location::gpu, access::mode::read).raw();

            computational.patchesState = patchesPd->getState(access::location::gpu, access::mode::read).raw();

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////

            computational.startType = storage.startType;
            computational.endType   = storage.endType;

            return computational;
        }

        //index_i is the current particle
        static inline __device__ EnergyForceTorque energyForceTorque(const int& index_i,const int& index_j,
                                                                     const ComputationalData& computational){

            Quat qs,qe;
            real3 lps,drse;

            int index_s = Helix_<potential>::setUpBond(index_i,index_j,
                                                       computational.patchesPos,
                                                       computational.patchesVector,
                                                       computational.patchesParentIndex,
                                                       computational.dir,
                                                       computational.box,
                                                       computational.startType,computational.endType,
                                                       qs,qe,lps,drse);

            const auto param = computational.paramPairIterator(index_i,index_j);

            const int state = computational.patchesState[index_s].z;

            //At this point we have the following information:
            // - q_s, q_e: the orientation of the patches
            // - lps: the vector from the center of the particle to the patch for the start particle
            // - drse: the vector from the patch of the start particle to the patch of the end particle
            // - state: the state of the start patch ( the state of the bond )

            ///////////////////////////////////////////////////

            real Eb;

            typename potential::params helixParams;

            tensor3 R_H;

            if(state == int(0)){

                Eb = param.Eb0;

                helixParams = param.helixParams0;

                R_H  = param.R_H0;

            } else {

                Eb = param.Eb1;

                helixParams = param.helixParams1;

                R_H  = param.R_H1;
            }

            if (index_s == index_j){
                R_H = R_H.transpose();
            }

            ///////////////////////////////////////////////////

            EnergyForceTorque eft;

            eft.energy = BasicPotentials::Helix::Dynamic::energy<potential>(drse,qs,qe,R_H,
                                                                            helixParams.dstParams,
                                                                            helixParams.thetaParams,
                                                                            helixParams.phiParams,
                                                                            Eb);

            eft.force  = make_real4(BasicPotentials::Helix::Dynamic::force<potential>
                                   (drse,qs,qe,R_H,
                                    helixParams.dstParams,
                                    helixParams.thetaParams,
                                    helixParams.phiParams,
                                    Eb),0.0);

            eft.torque = make_real4(BasicPotentials::Helix::Dynamic::torque<potential>
                                   (drse,lps,qs,qe,R_H,
                                    helixParams.dstParams,
                                    helixParams.thetaParams,
                                    helixParams.phiParams,
                                    Eb),0.0);

            return eft;
        }

        static inline __device__ StateTransitionProbability stateTransitionProbability(const int& index_i,const int& index_j,
                                                                                       const ComputationalData& computational){

            const real4 patch_i_pos = computational.patchesPos[index_i];
            const int tpy_i = int(patch_i_pos.w);

            StateTransitionProbability stp;

            if( tpy_i == computational.startType ){

                const int state = computational.patchesState[index_i].z;

                const auto param = computational.paramPairIterator(index_i,index_j);

                if(state == int(0)){
                    stp.tentativeState   = computational.patchesState[index_i];
                    stp.tentativeState.z = int(1);
                    stp.transitionProbability = param.prob_0_to_1;
                } else {
                    stp.tentativeState   = computational.patchesState[index_i];
                    stp.tentativeState.z = int(0);
                    stp.transitionProbability = param.prob_1_to_0;
                }

            } else {
                stp.tentativeState        = computational.patchesState[index_i];
                stp.transitionProbability = real(0.0);
            }

            return stp;
        }

    };

    using HelixExponential        = NonBondedDynamicallyBondedPatchyParticles_<Helix_<typename BasicPotentials::Helix::exponential>>;
    using HelixExponential2States = NonBondedDynamicallyBondedPatchyParticlesState_<Helix2States_<typename BasicPotentials::Helix::exponential>>;

    using HelixCosine        = NonBondedDynamicallyBondedPatchyParticles_<Helix_<typename BasicPotentials::Helix::cosine>>;
    using HelixCosine2States = NonBondedDynamicallyBondedPatchyParticlesState_<Helix2States_<typename BasicPotentials::Helix::cosine>>;

}}}}

REGISTER_NONBONDED_PATCHES_INTERACTOR(
    NonBondedPatches,HelixExponential,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBondedPatches::HelixExponential>
)

REGISTER_NONBONDED_PATCHES_INTERACTOR(
    NonBondedPatches,HelixExponential2States,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBondedPatches::HelixExponential2States>
)

REGISTER_NONBONDED_PATCHES_INTERACTOR(
    NonBondedPatches,HelixCosine,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBondedPatches::HelixCosine>
)

REGISTER_NONBONDED_PATCHES_INTERACTOR(
    NonBondedPatches,HelixCosine2States,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBondedPatches::HelixCosine2States>
)
