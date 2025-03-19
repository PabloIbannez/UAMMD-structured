#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace AFMimage{

template <class AFMimageType>
class AFMimage_{

    public:

        ///////////////////////////

        //Computational data
        using ComputationalData = typename AFMimageType::ComputationalData;

        //Potential parameters
        using StorageData = typename AFMimageType::StorageData;

        //Expose auxiliary functions
        static inline __device__ real OperationOverGenerated(const real& generatedPixelValue){
            return AFMimageType::OperationOverGenerated(generatedPixelValue);
        }

        static inline __device__ real OperationOverReferenceAndGenerated(const real& referencePixelValue,
                                                                         const real& generatedPixelValue){
            return AFMimageType::OperationOverReferenceAndGenerated(referencePixelValue,generatedPixelValue);
        }

        static inline __device__ real3 OperationOverReferenceAndGeneratedDerivative(const real& referencePixelValue,
                                                                                    const real3& generatedPixelValueDerivative){
            return AFMimageType::OperationOverReferenceAndGeneratedDerivative(referencePixelValue,generatedPixelValueDerivative);
        }

        static inline __device__ real3 OperationOverGeneratedAndGeneratedDerivative(const real& generatedPixelValue,
                                                                                    const real3& generatedPixelValueDerivative){
            return AFMimageType::OperationOverGeneratedAndGeneratedDerivative(generatedPixelValue,generatedPixelValueDerivative);
        }

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return AFMimageType::getComputationalData(this->gd,this->pg,
                                                      storage,comp,st);
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        AFMimage_(std::shared_ptr<GlobalData>    gd,
                  std::shared_ptr<ParticleGroup> pg,
                  DataEntry& data):gd(gd),
                                   pg(pg),pd(getExtendedParticleData(pg)){

            storage = AFMimageType::getStorageData(gd,pg,data);
        }

        //Transversers, energy, force

        struct EnergyTransverser{

            real*  energy;

            using resultType=real;

            EnergyTransverser(real* energy):energy(energy){}

            inline __device__ resultType zero(){
                return real(0.0);
            }

            inline __device__ void accumulate(resultType& total,const resultType current){
                total += current;
            }

            inline __device__ resultType compute(const int currentParticleIndex,
                                                 const ComputationalData& computational,
                                                 const real& opGen,const real& opRefGen,
                                                 const real3& opRefGenDer,const real3& opGenGenDer){

                return AFMimageType::energy(currentParticleIndex,computational,
                                            opGen,opRefGen,opRefGenDer,opGenGenDer);
            }

            inline __device__ void       set(const int& index_i,resultType& quantity){
                energy[index_i] += quantity;
            }

            inline __device__ resultType get(const int& index_i){
                return energy[index_i];
            }

        };

        struct ForceTransverser{

            real4* force;

            using resultType=real4;

            ForceTransverser(real4* force):force(force){}

            inline __device__ resultType zero(){
                return make_real4(0.0);
            }

            inline __device__ void accumulate(resultType& total,const resultType current){
                total += current;
            }

            inline __device__ resultType compute(const int currentParticleIndex,
                                                 const ComputationalData& computational,
                                                 const real& opGen,const real& opRefGen,
                                                 const real3& opRefGenDer,const real3& opGenGenDer){

                return make_real4(AFMimageType::force(currentParticleIndex,computational,
                                                      opGen,opRefGen,opRefGenDer,opGenGenDer),0.0);
            }

            inline __device__ void       set(const int& index_i,resultType& quantity){
                force[index_i] += quantity;
            }

            inline __device__ resultType get(const int& index_i){
                return force[index_i];
            }

        };

        EnergyTransverser getEnergyTransverser(){

            real* energy = this->pd->getEnergy(access::location::gpu,
                                              access::mode::readwrite).raw();
            return EnergyTransverser(energy);
        }

        ForceTransverser getForceTransverser(){

            real4* force = this->pd->getForce(access::location::gpu,
                                               access::mode::readwrite).raw();
            return ForceTransverser(force);
        }
};

}}}}
