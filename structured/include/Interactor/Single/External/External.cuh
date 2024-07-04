#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace External{

template <class ExternalType_>
struct EnergyTransverser_{

    real*  energy;

    using ExternalType = ExternalType_;
    using resultType   = real;

    EnergyTransverser_(real*  energy):energy(energy){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const typename ExternalType::ComputationalData& computational){
        return ExternalType::energy(index_i, computational);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        energy[index_i] += quantity;
    }
};

template <class ExternalType_>
struct ForceTransverser_{

    real4*  force;

    using ExternalType = ExternalType_;
    using resultType   = real4;

    ForceTransverser_(real4*  force):force(force){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const typename ExternalType::ComputationalData& computational){
        return make_real4(ExternalType::force(index_i, computational),0.0);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i] += quantity;
    }
};

template <class ExternalType_>
struct ForceTorqueTransverser_{

    real4*  force;
    real4*  torque;

    using ExternalType = ExternalType_;
    using resultType   = ForceTorque;

    ForceTorqueTransverser_(real4*  force,real4*  torque):force(force),torque(torque){}

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
                                         const typename ExternalType::ComputationalData& computational){
        return ExternalType::forceTorque(index_i, computational);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i]  += quantity.force;
        torque[index_i] += quantity.torque;
    }
};

  template <class ExternalType_>
  struct LambdaTransverser_{

    real* lambdaDerivative;

    using ExternalType  = ExternalType_;
    using resultType = real;

    LambdaTransverser_(real* lambdaDerivative):lambdaDerivative(lambdaDerivative){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const typename ExternalType::ComputationalData& computational){
      return ExternalType::lambdaDerivative(index_i, computational);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
      lambdaDerivative[index_i] += quantity;
    }
  };



  template <class ExternalType_>
  struct HessianTransverser_{

    tensor3*  hessian;
    const int*  id;
    const int*  selectedId;
    const int*  id2index;

    using ExternalType = ExternalType_;
    using resultType   = tensor3;

    HessianTransverser_(tensor3*  hessian, int* id, int* selectedId, int* id2index):hessian(hessian),id(id),
										    selectedId(selectedId),
										    id2index(id2index){}

    inline __device__ resultType zero(){return tensor3(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
      total  += current;
    }

    inline __device__ resultType compute(const int& index_i,
                                         const typename ExternalType::ComputationalData& computational){
      tensor3 H = tensor3(0.0);

      const int id_i  = id[index_i];
      const int selId = selectedId[index_i];

      if (selId == id_i){
	H = ExternalType::hessian(id_i, computational);
      }
      return H;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
      hessian[index_i]  += quantity;
    }
  };

template <class ExternalType_>
struct MagneticFieldTransverser_{

    real4* magneticField;

    using ExternalType = ExternalType_;
    using resultType    = real4;

    MagneticFieldTransverser_(real4* magneticField):magneticField(magneticField){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const typename ExternalType::ComputationalData& computational){
        return ExternalType::magneticField(index_i, computational);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        magneticField[index_i] += quantity;
    }
};

template <class ExternalType_>
struct ForceTorqueMagneticFieldTransverser_{

    real4*  force;
    real4*  torque;
    real4*  magneticField;

    using ExternalType = ExternalType_;
    using resultType    = ForceTorqueMagneticField;

    ForceTorqueMagneticFieldTransverser_(real4*  force,real4*  torque,real4*  magneticField)
      :force(force),torque(torque),magneticField(magneticField){}

    inline __device__ resultType zero(){
        resultType result;
        result.force        = make_real4(0.0);
        result.torque       = make_real4(0.0);
        result.magneticField= make_real4(0.0);
        return result;
    }

    inline __device__ void accumulate(resultType& total,const resultType current){
        total.force         += current.force;
        total.torque        += current.torque;
        total.magneticField += current.magneticField;
    }

    inline __device__ resultType compute(const int& index_i,
                                         const typename ExternalType::ComputationalData& computational){

        ForceTorque forceTorque = ExternalType::forceTorque(index_i, computational);
        real4    magneticField  = ExternalType::magneticField(index_i, computational);

        resultType result;

        result.force         = forceTorque.force;
        result.torque        = forceTorque.torque;
        result.magneticField = magneticField;

        return result;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i]         += quantity.force;
        torque[index_i]        += quantity.torque;
        magneticField[index_i] += quantity.magneticField;
    }
};

template<class ExternalType_>
class ExternalBase_{

    public:

        using ExternalType = ExternalType_;

        //Computational data
        using ComputationalData = typename ExternalType::ComputationalData;

        //Potential parameters
        using StorageData       = typename ExternalType::StorageData;

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return ExternalType::getComputationalData(this->gd,this->pg,
                                                      storage,comp,st);
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        ExternalBase_(std::shared_ptr<GlobalData>    gd,
                       std::shared_ptr<ParticleGroup> pg,
                       DataEntry& data):gd(gd),
                                        pg(pg),pd(getExtendedParticleData(pg)){

            storage = ExternalType::getStorageData(gd,pg,data);

        }

};

template<class ExternalType_>
class External_ : public ExternalBase_<ExternalType_>{

    public:

        using ExternalType = typename ExternalBase_<ExternalType_>::ExternalType;

        ///////////////////////////

        //Transversers

        using EnergyTransverser = EnergyTransverser_<ExternalType>;
        using ForceTransverser  = ForceTransverser_<ExternalType>;

        ///////////////////////////

        External_(std::shared_ptr<GlobalData>    gd,
                   std::shared_ptr<ParticleGroup> pg,
                   DataEntry& data):ExternalBase_<ExternalType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();

            return EnergyTransverser(energy);
        }

        ForceTransverser getForceTransverser(){

            real4*  force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();

            return ForceTransverser(force);
        }
  };

  template<class ExternalType_>
  class ExternalLambda_ : public External_<ExternalType_>{

  public:

    using ExternalType = typename External_<ExternalType_>::ExternalType;

    ///////////////////////////

    //Transversers

    using LambdaTransverser = LambdaTransverser_<ExternalType>;

    ///////////////////////////

    ExternalLambda_(std::shared_ptr<GlobalData>    gd,
		    std::shared_ptr<ParticleGroup> pg,
		    DataEntry& data):External_<ExternalType_>(gd,pg,data){}

    ///////////////////////////

    LambdaTransverser getLambdaTransverser(){

      real*  lambdaDerivative = this->pd->getLambdaDerivative(access::location::gpu, access::mode::readwrite).raw();

      return LambdaTransverser(lambdaDerivative);
    }
  };

  template<class ExternalType_>
  class ExternalTorque_ : public ExternalBase_<ExternalType_>{

    public:

        using ExternalType = typename ExternalBase_<ExternalType_>::ExternalType;

        ///////////////////////////

        //Transversers

        using EnergyTransverser = EnergyTransverser_<ExternalType>;
        using ForceTransverser  = ForceTorqueTransverser_<ExternalType>;

        ///////////////////////////

        ExternalTorque_(std::shared_ptr<GlobalData>    gd,
			std::shared_ptr<ParticleGroup> pg,
			DataEntry& data):ExternalBase_<ExternalType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();

            return EnergyTransverser(energy);
        }

        ForceTransverser getForceTransverser(){

            real4*  force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();

	    real4*  torque      = this->pd->getTorque(access::location::gpu, access::mode::readwrite).raw();
            return ForceTransverser(force, torque);
        }
};

  template<class ExternalType_>
  class ExternalHessian_ : public External_<ExternalType_>{

    public:

        using ExternalType = typename ExternalBase_<ExternalType_>::ExternalType;

        ///////////////////////////

        //Transversers

        using HessianTransverser = HessianTransverser_<ExternalType>;

        ///////////////////////////

        ExternalHessian_(std::shared_ptr<GlobalData>    gd,
			 std::shared_ptr<ParticleGroup> pg,
			 DataEntry& data):External_<ExternalType_>(gd,pg,data){}

        ///////////////////////////

        HessianTransverser getHessianTransverser(){

            tensor3*  hessian     = this->pd->getHessian(access::location::gpu, access::mode::readwrite).raw();
	    const int* id         = this->pd->getId(access::location::gpu, access::mode::read).raw();
	    const int* selectedId = this->pd->getSelectedId(access::location::gpu, access::mode::read).raw();
	    const int* id2index   = this->pd->getIdOrderedIndices(access::location::gpu);

	    return HessianTransverser(hessian, id, selectedId, id2index);
        }
};

template<class ExternalType_>
class ExternalForceTorqueMagneticField_ : public ExternalBase_<ExternalType_>{

    public:

        using ExternalType = typename ExternalBase_<ExternalType_>::ExternalType;

        ///////////////////////////

        //Transversers

        using EnergyTransverser             = EnergyTransverser_<ExternalType>;
        using ForceTransverser              = ForceTorqueTransverser_<ExternalType>;
        using ForceTorqueMagneticFieldTransverser = ForceTorqueMagneticFieldTransverser_<ExternalType>;
        using MagneticFieldTransverser      = MagneticFieldTransverser_<ExternalType>;

        ///////////////////////////

        ExternalForceTorqueMagneticField_(std::shared_ptr<GlobalData>    gd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          DataEntry& data):ExternalBase_<ExternalType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();

            return EnergyTransverser(energy);
        }

        ForceTransverser getForceTransverser(){

            real4*  force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            real4*  torque = this->pd->getTorque(access::location::gpu, access::mode::readwrite).raw();

            return ForceTransverser(force,torque);
        }

        ForceTorqueMagneticFieldTransverser getForceTorqueMagneticFieldTransverser(){

            real4*  force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            real4*  torque = this->pd->getTorque(access::location::gpu, access::mode::readwrite).raw();
            real4*  magneticField = this->pd->getMagneticField(access::location::gpu, access::mode::readwrite).raw();

            return ForceTorqueMagneticFieldTransverser(force,torque,magneticField);
        }

        MagneticFieldTransverser getMagneticFieldTransverser(){

            real4*  magneticField = this->pd->getMagneticField(access::location::gpu, access::mode::readwrite).raw();

            return MagneticFieldTransverser(magneticField);
        }
};

}}}}
