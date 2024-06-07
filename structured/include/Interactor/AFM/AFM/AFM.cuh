#ifndef __AFM__
#define __AFM__

namespace uammd{
namespace structured{
namespace Potentials{
namespace AFM{

template <class AFMType>
class AFM_{

    public:

        ///////////////////////////

        //Computational data
        using ComputationalData = typename AFMType::ComputationalData;

        //Potential parameters
        using StorageData = typename AFMType::StorageData;

        //AFM parameters
        using AFMParameters      = typename AFMType::AFMParameters;

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return AFMType::getComputationalData(this->gd,this->pg,
                                                 storage,comp,st);
        }

        template<typename T>
        AFMParameters processAFMParameters(std::map<std::string,T>& afmParametersMap){
            return AFMType::processAFMParameters(afmParametersMap);
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        AFM_(std::shared_ptr<GlobalData>    gd,
             std::shared_ptr<ParticleGroup> pg,
             DataEntry& data):gd(gd),
                              pg(pg),pd(getExtendedParticleData(pg)){

            storage = AFMType::getStorageData(gd,pg,data);
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
                                                 const int tipIndex,
                                                 const ComputationalData& computational,
                                                 const AFMParameters&      afmParam){

                return AFMType::energy(currentParticleIndex,tipIndex,computational,afmParam);
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
                                                 const int tipIndex,
                                                 const ComputationalData& computational,
                                                 const AFMParameters&      afmParam){

                return make_real4(AFMType::force(currentParticleIndex,tipIndex,computational,afmParam),0.0);
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

#endif
