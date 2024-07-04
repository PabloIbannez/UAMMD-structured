#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

template <class BondType_>
struct EnergyTransverser_{

    real*  energy;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = real;

    EnergyTransverser_(real*  energy,const int* id2index):energy(energy),id2index(id2index){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&   bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        return BondType::energy(i,j,currentParticleIndex,computational,bondParam)/real(2.0); //Each bond is counted twice
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        energy[index_i] += quantity;
    }
};

template <class BondType_>
struct ForceTransverser_{

    real4*  force;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = real4;

    ForceTransverser_(real4*  force,const int* id2index):force(force),id2index(id2index){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        return make_real4(BondType::force(i,j,currentParticleIndex,computational,bondParam),0.0);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i] += quantity;
    }
};

template <class BondType_>
struct LambdaTransverser_{

    real* lambdaDerivative;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = real;

    LambdaTransverser_(real* lambdaDerivative,const int* id2index):lambdaDerivative(lambdaDerivative),id2index(id2index){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        return BondType::lambdaDerivative(i,j,currentParticleIndex,computational,bondParam)/real(2.0); //Each bond is counted twice
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        lambdaDerivative[index_i] += quantity;
    }
};

template <class BondType_>
struct StressTransverser_{

    tensor3*  stress;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = tensor3;

    StressTransverser_(tensor3* stress,const int* id2index):stress(stress),id2index(id2index){}

    inline __device__ resultType zero(){return tensor3(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
        total+=current;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];

        real3 fij = BondType::force(i,j,currentParticleIndex,computational,bondParam); //This is the force over currentParticleIndex
                                                                                       //due to the other particle in the bond.

        real3 rij = computational.box.apply_pbc(make_real3(computational.pos[j]) - make_real3(computational.pos[i]));
        if        (currentParticleIndex == i){
        } else if (currentParticleIndex == j){
            rij = -rij;
        }

        return computeStress(rij,fij);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        stress[index_i] += quantity;
    }
};

template <class BondType_>
struct HessianTransverser_{

    tensor3*    hessian;
    const int*  id;
    const int*  selectedId;
    const int*  id2index;

    using BondType   = BondType_;
    using resultType = tensor3;

    HessianTransverser_(tensor3* hessian,
                        const int* id,
                        const int* selectedId,
                        const int* id2index):hessian(hessian),
                                             id(id),
                                             selectedId(selectedId),
                                             id2index(id2index){}

    inline __device__ resultType zero(){return tensor3(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
        total+=current;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){

        // 11 12 13
        // 21 22 23
        // 31 32 33

        //We compute a column of the hessian
        //The element hessian[i] will store the block (i,selectedId)
        //We first derevite i and then selectedId

        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];

        const int currentId = id[currentParticleIndex];
        const int selId     = selectedId[currentParticleIndex];

        tensor3 H = tensor3(0.0);

        if        (currentParticleIndex == i){
            if (bondParam.id_j == selId){
                H = BondType::hessian(i,j,currentParticleIndex,computational,bondParam);
            }
        } else if (currentParticleIndex == j){
            if (bondParam.id_i == selId){
	      H = BondType::hessian(i,j,currentParticleIndex,computational,bondParam);
            }
        }

        //Self term
        if(currentId == selId){
	  H = real(-1.0)*BondType::hessian(i,j,currentParticleIndex,computational,bondParam);
        }

        return H;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        hessian[index_i] += quantity;
    }
};

  template <class BondType_>
  struct PairwiseForceTransverser_{

    real4*    pairwiseForce;
    const int*  id;
    const int*  selectedId;
    const int*  id2index;

    using BondType   = BondType_;
    using resultType = real4;

    PairwiseForceTransverser_(real4* pairwiseForce,
			       const int* id,
			       const int* selectedId,
			       const int* id2index):pairwiseForce(pairwiseForce),
						    id(id),
						    selectedId(selectedId),
						    id2index(id2index){}

    inline __device__ resultType zero(){return real4();}

    inline __device__ void accumulate(resultType& total,const resultType current){
      total+=current;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){

      const int i = id2index[bondParam.id_i];
      const int j = id2index[bondParam.id_j];

      const int currentId = id[currentParticleIndex];
      const int selId     = selectedId[currentParticleIndex];

      real4 f = real4();
      if ((currentParticleIndex == i and bondParam.id_j == selId) or (currentParticleIndex == j and bondParam.id_i == selId)){
	f = -make_real4(BondType::force(i,j,currentParticleIndex,computational,bondParam),0.0); //Force that particle currentId makes over particle selId
      }
      return f;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
      pairwiseForce[index_i] += quantity;
    }
  };

  template <class BondType_>
  struct ForceTorqueTransverser_{

    real4*  force;
    real4*  torque;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = ForceTorque;

    ForceTorqueTransverser_(real4*  force,
                            real4*  torque,
                            const int* id2index):force(force),
                                                 torque(torque),
                                                 id2index(id2index){}

    inline __device__ resultType zero(){
        resultType zero;
        zero.force  = make_real4(0.0);
        zero.torque = make_real4(0.0);
        return zero;
    }

    inline __device__ void accumulate(resultType& total,const resultType current){
        total.force  += current.force;
        total.torque += current.torque;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        return BondType::forceTorque(i,j,currentParticleIndex,computational,bondParam);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i]  += quantity.force;
        torque[index_i] += quantity.torque;
    }
};

template<class BondType_>
class Bond2Base_ {

    public:

        ///////////////////////////

        //Number of particles in the bond type

        static constexpr int nPart = 2;

        std::vector<std::string> getParticleBondLabels(){
            std::vector<std::string> labels = {"id_i","id_j"};
            return labels;
        }

        ///////////////////////////

        struct BondType : public BondType_{
            //Bond parameters
            struct BondParameters : public BondType_::BondParameters {
                int id_i;
                int id_j;
            };
        };

        ///////////////////////////

        //Computational data
        using ComputationalData = typename BondType::ComputationalData;

        //Potential parameters
        using StorageData       = typename BondType::StorageData;

        //Bond parameters
        using BondParameters    = typename BondType::BondParameters;

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return BondType::getComputationalData(this->gd,
                                                  this->pg,storage,comp,st);
        }

        template<typename T>
        BondParameters processBondParameters(std::map<std::string,T>& bondParametersMap){
            BondParameters param;

            static_cast<typename BondType_::BondParameters&>(param) = BondType_::processBondParameters(this->gd,bondParametersMap);
            param.id_i = bondParametersMap.at("id_i");
            param.id_j = bondParametersMap.at("id_j");

            return param;
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        Bond2Base_(std::shared_ptr<GlobalData>    gd,
                   std::shared_ptr<ParticleGroup> pg,
                   DataEntry& data):gd(gd),
                                    pg(pg),pd(getExtendedParticleData(pg)){

            storage = BondType::getStorageData(gd,pg,data);
        }

};

template<class BondType_>
class Bond2_ : public Bond2Base_<BondType_>{

    public:

        using BondType = typename Bond2_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser        = EnergyTransverser_<BondType>;
        using ForceTransverser         = ForceTransverser_<BondType>;
        using StressTransverser        = StressTransverser_<BondType>;
        using PairwiseForceTransverser = PairwiseForceTransverser_<BondType>;

        ///////////////////////////

        Bond2_(std::shared_ptr<GlobalData>    gd,
               std::shared_ptr<ParticleGroup> pg,
               DataEntry& data):Bond2Base_<BondType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy       = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(energy,id2index);
        }

        ForceTransverser getForceTransverser(){

            real4*  force       = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(force,id2index);
        }

        StressTransverser getStressTransverser(){

            tensor3*  stress    = this->pd->getStress(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return StressTransverser(stress,id2index);
        }

        PairwiseForceTransverser getPairwiseForceTransverser(){

            real4*  pforce        = this->pd->getPairwiseForce(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index   = this->pd->getIdOrderedIndices(access::location::gpu);
	    const int* id         = this->pd->getId(access::location::gpu, access::mode::read).raw();
            const int* selectedId = this->pd->getSelectedId(access::location::gpu, access::mode::read).raw();


            return PairwiseForceTransverser(pforce,
					    id,
					    selectedId,id2index);
        }
};

template<class BondType_>
class Bond2Lambda_ : public Bond2_<BondType_>{

    public:

        using BondType = typename Bond2_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using LambdaTransverser = LambdaTransverser_<BondType>;

        Bond2Lambda_(std::shared_ptr<GlobalData>    gd,
                     std::shared_ptr<ParticleGroup> pg,
                     DataEntry& data):Bond2_<BondType_>(gd,pg,data){}

        ///////////////////////////

        LambdaTransverser getLambdaTransverser(){

            real*  lambdaDerivative = this->pd->getLambdaDerivative(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index     = this->pd->getIdOrderedIndices(access::location::gpu);

            return LambdaTransverser(lambdaDerivative,id2index);
        }
};

template<class BondType_>
class Bond2Hessian_ : public Bond2_<BondType_> {

    public:

        using BondType = typename Bond2_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using HessianTransverser = HessianTransverser_<BondType>;

    public:

        Bond2Hessian_(std::shared_ptr<GlobalData>    gd,
                      std::shared_ptr<ParticleGroup> pg,
                      DataEntry& data):Bond2_<BondType_>(gd,pg,data){}

        ///////////////////////////

        HessianTransverser getHessianTransverser(){

            tensor3*  hessian     = this->pd->getHessian(access::location::gpu, access::mode::readwrite).raw();
            const int* id         = this->pd->getId(access::location::gpu, access::mode::read).raw();
            const int* selectedId = this->pd->getSelectedId(access::location::gpu, access::mode::read).raw();
            const int* id2index   = this->pd->getIdOrderedIndices(access::location::gpu);

            return HessianTransverser(hessian,
                                      id,
                                      selectedId,id2index);
        }
};

template<class BondType_>
class Bond2Torque_ : public Bond2Base_<BondType_> {

    public:

        using BondType = typename Bond2_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = EnergyTransverser_<BondType>;
        using ForceTransverser  = ForceTorqueTransverser_<BondType>;

    public:

        Bond2Torque_(std::shared_ptr<GlobalData>    gd,
                     std::shared_ptr<ParticleGroup> pg,
                     DataEntry& data):Bond2Base_<BondType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy       = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(energy,id2index);
        }

        ForceTransverser getForceTransverser(){

            real4*  force       = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            real4*  torque      = this->pd->getTorque(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(force,torque,id2index);
        }
};

}}}}

