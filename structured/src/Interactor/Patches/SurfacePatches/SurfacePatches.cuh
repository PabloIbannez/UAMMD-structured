#pragma once

#include "UAMMDstructuredBase.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace SurfacePatches{

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
struct ForceTorqueTransverser_{

    real4*  force;
    real4*  torque;

    using SurfaceType = SurfaceType_;
    using resultType    = ForceTorque;

    ForceTorqueTransverser_(real4*  force,real4*  torque):force(force),torque(torque){}

    inline __device__ resultType zero(){
        ForceTorque result;
        result.force = make_real4(0.0);
        result.torque = make_real4(0.0);
        return result;
    }


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


template<class SurfaceType_>
class SurfacePatchyParticles_{

    public:

        using SurfaceType = SurfaceType_;

        //Transversers

        using EnergyTransverser = EnergyTransverser_<SurfaceType>;
        using ForceTransverser  = ForceTorqueTransverser_<SurfaceType>;

        //Computational data
        using ComputationalData = typename SurfaceType::ComputationalData;

        //Potential parameters
        using StorageData       = typename SurfaceType::StorageData;

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return SurfaceType::getComputationalData(this->gd, this->pg,
                                                     this->patchesGd,this->patchesPg,
                                                     storage,comp,st);
        }

    private:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<ExtendedParticleData> pd;

        std::shared_ptr<GlobalData>    patchesGd;
        std::shared_ptr<ParticleGroup> patchesPg;
        std::shared_ptr<ExtendedParticleData> patchesPd;

        StorageData storage;

    public:

        SurfacePatchyParticles_(std::shared_ptr<GlobalData>    gd,
                                std::shared_ptr<ParticleGroup> pg,
                                std::shared_ptr<GlobalData>    patchesGd,
                                std::shared_ptr<ParticleGroup> patchesPg,
                                DataEntry& data):gd(gd),              pg(pg),              pd(getExtendedParticleData(pg)),
                                                 patchesGd(patchesGd),patchesPg(patchesPg),patchesPd(getExtendedParticleData(patchesPg)){

            storage = SurfaceType::getStorageData(gd,pg,patchesGd,patchesPg,data);
        }

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  patchesEnergy = patchesPd->getEnergy(access::location::gpu, access::mode::readwrite).raw();

            return EnergyTransverser(patchesEnergy);
        }

        ForceTransverser getForceTransverser(){

            real4* patchesForce  = patchesPd->getForce(access::location::gpu, access::mode::readwrite).raw();
            real4* patchesTorque = patchesPd->getTorque(access::location::gpu, access::mode::readwrite).raw();

            return ForceTransverser(patchesForce,patchesTorque);
        }
};

}}}}


