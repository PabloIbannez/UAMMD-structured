#pragma once

#include "UAMMDstructuredBase.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Surface{

template <class SurfaceType_>
struct EnergyTransverser_{

    real*  energy;

    using SurfaceType = SurfaceType_;
    using resultType   = real;

    EnergyTransverser_(real*  energy):energy(energy){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const typename SurfaceType::ComputationalData& computational){
        return SurfaceType::energy(index_i, computational);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        energy[index_i] += quantity;
    }
};

template <class SurfaceType_>
struct ForceTransverser_{

    real4*  force;

    using SurfaceType = SurfaceType_;
    using resultType   = real4;

    ForceTransverser_(real4*  force):force(force){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const typename SurfaceType::ComputationalData& computational){
        return make_real4(SurfaceType::force(index_i, computational),0.0);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i] += quantity;
    }
};

template <class SurfaceType_>
struct ForceTorqueTransverser_{

    real4*  force;
    real4*  torque;

    using SurfaceType = SurfaceType_;
    using resultType    = ForceTorque;

    ForceTorqueTransverser_(real4*  force,real4*  torque):force(force),torque(torque){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
        total.force  += current.force;
        total.torque += current.torque;
    }

    inline __device__ resultType compute(const int& index_i,
                                         const typename SurfaceType::ComputationalData& computational){
        return SurfaceType::forceTorque(index_i, computational);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i]  += quantity.force;
        torque[index_i] += quantity.torque;
    }
};

  template <class SurfaceType_>
  struct HessianTransverser_{

    tensor3*  hessian;
    const int*  id;
    const int*  selectedId;
    const int*  id2index;

    using SurfaceType = SurfaceType_;
    using resultType    = tensor3;

    HessianTransverser_(tensor3*  hessian, int* id, int* selectedId, int* id2index):hessian(hessian),id(id),
										    selectedId(selectedId),
										    id2index(id2index){}

    inline __device__ resultType zero(){return tensor3(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
      total  += current;
    }

    inline __device__ resultType compute(const int& index_i,
                                         const typename SurfaceType::ComputationalData& computational){
      tensor3 H = tensor3(0.0);

      const int id_i  = id[index_i];
      const int selId = selectedId[index_i];

      if (selId == id_i){
	H = SurfaceType::hessian(index_i, computational);
      }
      return H;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
      hessian[index_i]  += quantity;
    }
  };

template <class SurfaceType_>
struct MagneticFieldTransverser_{

    real4* magneticField;

    using SurfaceType = SurfaceType_;
    using resultType    = real4;

    MagneticFieldTransverser_(real4* magneticField):magneticField(magneticField){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int& index_i,
                                         const typename SurfaceType::ComputationalData& computational){
        return SurfaceType::magneticField(index_i, computational);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        magneticField[index_i] += quantity;
    }
};

template<class SurfaceType_>
class SurfaceBase_{

    public:

        using SurfaceType = SurfaceType_;

        //Computational data
        using ComputationalData = typename SurfaceType::ComputationalData;

        //Potential parameters
        using StorageData       = typename SurfaceType::StorageData;

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return SurfaceType::getComputationalData(this->gd,this->pg,
                                                     storage,comp,st);
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        SurfaceBase_(std::shared_ptr<GlobalData>    gd,
                     std::shared_ptr<ParticleGroup> pg,
                     DataEntry& data):gd(gd),
                                      pg(pg),pd(getExtendedParticleData(pg)){

            storage = SurfaceType::getStorageData(gd,pg,data);

        }

};

template<class SurfaceType_>
class Surface_ : public SurfaceBase_<SurfaceType_>{

    public:

        using SurfaceType = typename SurfaceBase_<SurfaceType_>::SurfaceType;

        ///////////////////////////

        //Transversers

        using EnergyTransverser = EnergyTransverser_<SurfaceType>;
        using ForceTransverser  = ForceTransverser_<SurfaceType>;

        ///////////////////////////

        Surface_(std::shared_ptr<GlobalData>    gd,
                 std::shared_ptr<ParticleGroup> pg,
                 DataEntry& data):SurfaceBase_<SurfaceType_>(gd,pg,data){}

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


  template<class SurfaceType_>
  class SurfaceHessian_ : public Surface_<SurfaceType_>{

    public:

        using SurfaceType = typename SurfaceBase_<SurfaceType_>::SurfaceType;

        ///////////////////////////

        //Transversers

        using HessianTransverser = HessianTransverser_<SurfaceType>;

        ///////////////////////////

        SurfaceHessian_(std::shared_ptr<GlobalData>    gd,
			 std::shared_ptr<ParticleGroup> pg,
			 DataEntry& data):Surface_<SurfaceType_>(gd,pg,data){}

        ///////////////////////////

        HessianTransverser getHessianTransverser(){

            tensor3*  hessian     = this->pd->getHessian(access::location::gpu, access::mode::readwrite).raw();
	    const int* id         = this->pd->getId(access::location::gpu, access::mode::read).raw();
	    const int* selectedId = this->pd->getSelectedId(access::location::gpu, access::mode::read).raw();
	    const int* id2index   = this->pd->getIdOrderedIndices(access::location::gpu);

	    return HessianTransverser(hessian, id, selectedId, id2index);
        }
};

template<class SurfaceType_>
class SurfaceTorque_ : public SurfaceBase_<SurfaceType_>{

    public:

        using SurfaceType = typename SurfaceBase_<SurfaceType_>::SurfaceType;

        ///////////////////////////

        //Transversers

        using EnergyTransverser = EnergyTransverser_<SurfaceType>;
        using ForceTransverser  = ForceTorqueTransverser_<SurfaceType>;

        ///////////////////////////

        SurfaceTorque_(std::shared_ptr<GlobalData>    gd,
                       std::shared_ptr<ParticleGroup> pg,
                       DataEntry& data):SurfaceBase_<SurfaceType_>(gd,pg,data){}

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
};

template<class SurfaceType_>
class SurfaceTorqueMagneticField_ : public SurfaceTorque_<SurfaceType_>{

    public:

        using SurfaceType = typename SurfaceTorque_<SurfaceType_>::SurfaceType;

        ///////////////////////////

        //Transversers

        using MagneticFieldTransverser = MagneticFieldTransverser_<SurfaceType>;

        ///////////////////////////

        SurfaceTorqueMagneticField_(std::shared_ptr<GlobalData>    gd,
                                    std::shared_ptr<ParticleGroup> pg,
                                    DataEntry& data):SurfaceTorque_<SurfaceType_>(gd,pg,data){}

        ///////////////////////////

        MagneticFieldTransverser getMagneticFieldTransverser(){

            real4*  magneticField = this->pd->getMagneticField(access::location::gpu, access::mode::readwrite).raw();

            return MagneticFieldTransverser(magneticField);
        }
};

}}}}

