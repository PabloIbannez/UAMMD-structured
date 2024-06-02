#ifndef __BOND1__
#define __BOND1__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond1{

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
        return BondType::energy(i,currentParticleIndex,computational,bondParam);
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
        return make_real4(BondType::force(i,currentParticleIndex,computational,bondParam),0.0);
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
        return BondType::lambdaDerivative(i,currentParticleIndex,computational,bondParam);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        lambdaDerivative[index_i] += quantity;
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

        const int currentId = id[currentParticleIndex];
        const int selId     = selectedId[currentParticleIndex];

        tensor3 H = tensor3(0.0);

        if (selId == currentId){
	  H = BondType::hessian(i,currentParticleIndex,computational,bondParam);
	}
	return H;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        hessian[index_i] += quantity;
    }
};


template <class BondType_>
struct ForceTorqueTransverser_{

    real4*  force;
    real4*  torque;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = ForceTorque;

    ForceTorqueTransverser_(real4*  force,real4*  torque,const int* id2index):force(force),torque(torque),id2index(id2index){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
        total.force  += current.force;
        total.torque += current.torque;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        return BondType::forceTorque(i,currentParticleIndex,computational,bondParam);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i]  += quantity.force;
        torque[index_i] += quantity.torque;
    }
};

template<class BondType_>
class Bond1Base_ {

    public:

        ///////////////////////////

        //Number of particles in the bond type

        static constexpr int nPart = 1;

        std::vector<std::string> getParticleBondLabels(){
            std::vector<std::string> labels = {"id_i"};
            return labels;
        }

        ///////////////////////////

        struct BondType : public BondType_{
            //Bond parameters
            struct BondParameters : public BondType_::BondParameters {
                int id_i;
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

        ComputationalData getComputationalData(const Computables& comp){
            return BondType::getComputationalData(this->gd,
                                                  this->pg,storage,comp);
        }

        template<typename T>
        BondParameters processBondParameters(std::map<std::string,T>& bondParametersMap){
            BondParameters param;

            static_cast<typename BondType_::BondParameters&>(param) = BondType_::processBondParameters(this->gd,bondParametersMap);
            param.id_i = bondParametersMap.at("id_i");

            return param;
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        Bond1Base_(std::shared_ptr<GlobalData>    gd,
                   std::shared_ptr<ParticleGroup> pg,
                   DataEntry& data):gd(gd),
                                    pg(pg),pd(getExtendedParticleData(pg)){

            storage = BondType::getStorageData(gd,pg,data);
        }

};

template<class BondType_>
class Bond1_ : public Bond1Base_<BondType_>{

    public:

        using BondType = typename Bond1_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = EnergyTransverser_<BondType>;
        using ForceTransverser  = ForceTransverser_<BondType>;

        ///////////////////////////

        Bond1_(std::shared_ptr<GlobalData>    gd,
               std::shared_ptr<ParticleGroup> pg,
               DataEntry& data):Bond1Base_<BondType_>(gd,pg,data){}

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
};

  template<class BondType_>
  class Bond1Hessian_ : public Bond1_<BondType_> {

    public:

        using BondType = typename Bond1_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using HessianTransverser = HessianTransverser_<BondType>;

        ///////////////////////////

        Bond1Hessian_(std::shared_ptr<GlobalData>    gd,
		      std::shared_ptr<ParticleGroup> pg,
		      DataEntry& data):Bond1_<BondType_>(gd,pg,data){}

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
class Bond1Lambda_ : public Bond1_<BondType_> {

    public:

        using BondType = typename Bond1_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using LambdaTransverser = LambdaTransverser_<BondType>;

        ///////////////////////////

        Bond1Lambda_(std::shared_ptr<GlobalData>    gd,
                     std::shared_ptr<ParticleGroup> pg,
                     DataEntry& data):Bond1_<BondType_>(gd,pg,data){}

        LambdaTransverser getLambdaTransverser(){

            real*  lambdaDerivative = this->pd->getLambdaDerivative(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index     = this->pd->getIdOrderedIndices(access::location::gpu);

            return LambdaTransverser(lambdaDerivative,id2index);
        }

};

template<class BondType_>
class Bond1Torque_ : public Bond1Base_<BondType_> {

    public:

        using BondType = typename Bond1_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = EnergyTransverser_<BondType>;
        using ForceTransverser  = ForceTorqueTransverser_<BondType>;

    public:

        Bond1Torque_(std::shared_ptr<GlobalData>    gd,
                     std::shared_ptr<ParticleGroup> pg,
                     DataEntry& data):Bond1Base_<BondType_>(gd,pg,data){}

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
};

}}}}

#endif
