#ifndef __NON_BONDED_PATCHES_POTENTIAL__
#define __NON_BONDED_PATCHES_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBondedPatches{

template <class NonBondedType_>
struct EnergyTransverser_{

    int*   patchesId;

    int4*  patchesState;
    real4* patchesStateInfo;

    real*  patchesEnergy;

    using NonBondedType = NonBondedType_;
    using resultType    = real;

    EnergyTransverser_(int*   patchesId,
                       int4*  patchesState,
                       real4* patchesStateInfo,
                       real*  patchesEnergy):patchesId(patchesId),
                                             patchesState(patchesState),
                                             patchesStateInfo(patchesStateInfo),
                                             patchesEnergy(patchesEnergy){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const int& index_j,
                                         const typename NonBondedType::ComputationalData& computational){

        EnergyForceTorque eFrcTrq;

        eFrcTrq.energy = real(0.0);
        eFrcTrq.force  = make_real4(0.0);
        eFrcTrq.torque = make_real4(0.0);

        bool bonded = (patchesState[index_i].x == patchesId[index_j] and
                       patchesState[index_j].x == patchesId[index_i]);

        if((patchesState[index_i].x == int(-1) and
            patchesState[index_j].x == int(-1) )
            or bonded)
        {
            eFrcTrq = NonBondedType::energyForceTorque(index_i,index_j,computational);

            if(eFrcTrq.energy < patchesStateInfo[index_i].x){
                patchesStateInfo[index_i].x = eFrcTrq.energy;
                patchesState[index_i].y     = patchesId[index_j];
            }
        }

        return eFrcTrq.energy/real(2.0); // We divide by 2 because we are computing the energy twice
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        patchesEnergy[index_i] += quantity;
    }
};

template <class NonBondedType_>
struct ForceTorqueTransverser_{

    int*   patchesId;

    int4*  patchesState;
    real4* patchesStateInfo;

    real4*  patchesForce;
    real4*  patchesTorque;

    using NonBondedType = NonBondedType_;
    using resultType    = ForceTorque;

    ForceTorqueTransverser_(int*   patchesId,
                            int4*  patchesState,
                            real4* patchesStateInfo,
                            real4* patchesForce,
                            real4* patchesTorque):patchesId(patchesId),
                                                  patchesState(patchesState),
                                                  patchesStateInfo(patchesStateInfo),
                                                  patchesForce(patchesForce),
                                                  patchesTorque(patchesTorque){}

    inline __device__ resultType zero(){
        resultType result;
        result.force  = make_real4(0.0);
        result.torque = make_real4(0.0);
        return result;
    }

    inline __device__ void accumulate(resultType& total,const resultType current){
        total.force  += current.force;
        total.torque += current.torque;
    }

    inline __device__ resultType compute(const int& index_i,
                                         const int& index_j,
                                         const typename NonBondedType::ComputationalData& computational){
        EnergyForceTorque eFrcTrq;

        eFrcTrq.energy = real(0.0);
        eFrcTrq.force  = make_real4(0.0);
        eFrcTrq.torque = make_real4(0.0);

        bool bonded = (patchesState[index_i].x == patchesId[index_j] and
                       patchesState[index_j].x == patchesId[index_i]);

        if((patchesState[index_i].x == int(-1) and
            patchesState[index_j].x == int(-1) )
            or bonded)
        {
            eFrcTrq = NonBondedType::energyForceTorque(index_i,index_j,computational);

            if(eFrcTrq.energy < patchesStateInfo[index_i].x){
                patchesStateInfo[index_i].x = eFrcTrq.energy;
                patchesState[index_i].y     = patchesId[index_j];
            }
        }

        resultType frcTrq;

        frcTrq.force  = eFrcTrq.force;
        frcTrq.torque = eFrcTrq.torque;

        return frcTrq;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        patchesForce[index_i]  += quantity.force;
        patchesTorque[index_i] += quantity.torque;
    }
};

template <class NonBondedType_>
struct TransitionProbabilityTransverser_{

    int*   patchesId;

    int4*  patchesState;
    real4* patchesStateInfo;

		int4* tentativeState;
		real* transProbability;

    using NonBondedType = NonBondedType_;
    using resultType    = StateTransitionProbability;

    TransitionProbabilityTransverser_(int*   patchesId,
                                      int4*  patchesState,
                                      real4* patchesStateInfo,
                                      int4*  tentativeState,
                                      real*  transProbability):patchesId(patchesId),
                                                               patchesState(patchesState),
                                                               patchesStateInfo(patchesStateInfo),
                                                               tentativeState(tentativeState),
                                                               transProbability(transProbability){}

    inline __device__ resultType zero(){
        resultType result;
        result.tentativeState        = {-1,-1,0,0};
        result.transitionProbability = real(0.0);
        return result;
    }

    inline __device__ void accumulate(resultType& total,const resultType current){
        total.transitionProbability = current.transitionProbability; //No need to accumulate, patch can only be bonded to one other patch
        total.tentativeState        = current.tentativeState;
    }

    inline __device__ resultType compute(const int& index_i,
                                         const int& index_j,
                                         const typename NonBondedType::ComputationalData& computational){
        resultType result;

        bool bonded = (patchesState[index_i].x == patchesId[index_j] and
                       patchesState[index_j].x == patchesId[index_i]);

        if(bonded){
            //If bonded, compute the transition probability
            result = NonBondedType::stateTransitionProbability(index_i,index_j,computational);
        } else {
            //Transition probability is 0.0
            result.transitionProbability = real(0.0);
            result.tentativeState        = patchesState[index_i];
        }

        return result;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        transProbability[index_i] = quantity.transitionProbability; //No need to accumulate, patch can only be bonded to one other patch
        tentativeState[index_i]   = quantity.tentativeState;
    }
};

template<class NonBondedType_>
class NonBondedDynamicallyBondedPatchyParticles_{

    public:

        using NonBondedType = NonBondedType_;

        ///////////////////////////

        //Computational data
        using ComputationalData = typename NonBondedType::ComputationalData;

        //Potential parameters
        using StorageData       = typename NonBondedType::StorageData;

        //Transversers

        using EnergyTransverser = EnergyTransverser_<NonBondedType>;
        using ForceTransverser  = ForceTorqueTransverser_<NonBondedType>;

        ComputationalData getComputationalData(){
            return NonBondedType::getComputationalData(this->gd, this->pg,
                                                        this->patchesGd,this->patchesPg,
                                                        storage);
        }

    private:

        std::shared_ptr<GlobalData>           gd;
        std::shared_ptr<ParticleGroup>        pg;
        std::shared_ptr<ExtendedParticleData> pd;

        std::shared_ptr<GlobalData>           patchesGd;
        std::shared_ptr<ParticleGroup>        patchesPg;
        std::shared_ptr<ExtendedParticleData> patchesPd;

        StorageData storage;

    public:

        NonBondedDynamicallyBondedPatchyParticles_(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   std::shared_ptr<GlobalData>    patchesGd,
                                                   std::shared_ptr<ParticleGroup> patchesPg,
                                                   DataEntry& data):gd(gd),              pg(pg),              pd(getExtendedParticleData(pg)),
                                                                    patchesGd(patchesGd),patchesPg(patchesPg),patchesPd(getExtendedParticleData(patchesPg)){

            storage = NonBondedType::getStorageData(gd,pg,patchesGd,patchesPg,data);
        }

        real getCutOff(){return storage.cutOff;}

        EnergyTransverser getEnergyTransverser(){

            int* patchesId    = patchesPd->getId(access::location::gpu, access::mode::read).raw();

            int4* patchesState      = patchesPd->getState(access::location::gpu, access::mode::read).raw();
            real4* patchesStateInfo = patchesPd->getStateInfo(access::location::gpu, access::mode::read).raw();

            real*  patchesEnergy = patchesPd->getEnergy(access::location::gpu, access::mode::readwrite).raw();

            return EnergyTransverser(patchesId,
                                     patchesState,
                                     patchesStateInfo,
                                     patchesEnergy);
        }

        ForceTransverser getForceTransverser(){

            int* patchesId    = patchesPd->getId(access::location::gpu, access::mode::read).raw();

            int4* patchesState      = patchesPd->getState(access::location::gpu, access::mode::read).raw();
            real4* patchesStateInfo = patchesPd->getStateInfo(access::location::gpu, access::mode::read).raw();

            real4* patchesForce  = patchesPd->getForce(access::location::gpu, access::mode::readwrite).raw();
            real4* patchesTorque = patchesPd->getTorque(access::location::gpu, access::mode::readwrite).raw();

            return ForceTransverser(patchesId,
                                    patchesState,
                                    patchesStateInfo,
                                    patchesForce,
                                    patchesTorque);
        }
};

template<class NonBondedType_>
class NonBondedDynamicallyBondedPatchyParticlesState_ :
    public NonBondedDynamicallyBondedPatchyParticles_<NonBondedType_>{

        public:

            using NonBondedType = typename NonBondedDynamicallyBondedPatchyParticles_<NonBondedType_>::NonBondedType;

            ///////////////////////////

            //Transversers

            using TransitionProbabilityTransverser = TransitionProbabilityTransverser_<NonBondedType>;

            ///////////////////////////

            NonBondedDynamicallyBondedPatchyParticlesState_(std::shared_ptr<GlobalData>    gd,
                                                            std::shared_ptr<ParticleGroup> pg,
                                                            std::shared_ptr<GlobalData>    patchesGd,
                                                            std::shared_ptr<ParticleGroup> patchesPg,
                                                            DataEntry& data):NonBondedDynamicallyBondedPatchyParticles_<NonBondedType_>(gd,pg,patchesGd,patchesPg,data){}

            TransitionProbabilityTransverser getTransitionProbabilityTransverser(){

                int* patchesId    = this->patchesPd->getId(access::location::gpu, access::mode::read).raw();

                int4* patchesState      = this->patchesPd->getState(access::location::gpu, access::mode::read).raw();
                real4* patchesStateInfo = this->patchesPd->getStateInfo(access::location::gpu, access::mode::read).raw();

                real*  transProbability = this->patchesPd->getTransProbability(access::location::gpu, access::mode::readwrite).raw();
                int4*  tentativeState   = this->patchesPd->getTentativeState(access::location::gpu, access::mode::readwrite).raw();

                return TransitionProbabilityTransverser(patchesId,
                                                         patchesState,
                                                         patchesStateInfo,
                                                         transProbability,
                                                         tentativeState);
            }
};

}}}}

#endif
