#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace Set1{

template <class PotentialType_,int nSet, int nReductions>
struct CenterOfMassEnergyTransverser_{

    real*  energy;

    using PotentialType = PotentialType_;

    using reduceType   = real4;
    using propertyType = real;
    using spreadType   = real;

    using ComputationalData = typename PotentialType::ComputationalData;
    using SetParameters     = typename PotentialType::SetParameters;

    CenterOfMassEnergyTransverser_(real*  energy):energy(energy){}

    inline __device__ reduceType zero(const ComputationalData& computational,
                                      const SetParameters&    setParam,
                                      const int& reductionStep,
                                      const int& setIndex,
                                      const reduceType reduced[nReductions][nSet]){
        return make_real4(0.0);
    }

    inline __device__ reduceType transform(const int& index ,
                                           const ComputationalData& computational,
                                           const SetParameters&    setParam,
                                           const int& reductionStep,
                                           const int& setIndex,
                                           const reduceType reduced[nReductions][nSet]){

        const real m = computational.mass[index];
        return make_real4(m*make_real3(computational.pos[index]),m);
    }

    static inline __device__ reduceType reduce(const reduceType& reduced1,
                                               const reduceType& reduced2){
        return reduced1+reduced2;
    }

    //After reduction

    inline __device__ propertyType compute(const ComputationalData& computational,
                                           const SetParameters&    setParam,
                                           const reduceType reduced[nReductions][nSet]){

        const real  totalMass = reduced[0][0].w;
        const real3 com = computational.box.apply_pbc(make_real3(reduced[0][0])/totalMass);

        return PotentialType::energy(com,computational,setParam);
    }

    inline __device__ spreadType spread(const propertyType& quantity,
                                        const int& index,
                                        const ComputationalData& computational,
                                        const SetParameters&    setParam,
                                        const int& setIndex,
                                        const reduceType reduced[nReductions][nSet]){

        const real totalMass = reduced[0][setIndex].w;

        return (computational.mass[index]/totalMass)*quantity;
    }

    inline __device__ void set(const int& index,
                               const spreadType& quantity,
                               const ComputationalData& computational,
                               const SetParameters&    setParam,
                               const int& setIndex,
                               const reduceType reduced[nReductions][nSet]){

        energy[index] += quantity;
    }
};

template <class PotentialType_,int nSet, int nReductions>
struct CenterOfMassForceTransverser_{

    real4*  force;

    using PotentialType = PotentialType_;

    using reduceType   = real4;
    using propertyType = real4;
    using spreadType   = real4;

    using ComputationalData = typename PotentialType::ComputationalData;
    using SetParameters     = typename PotentialType::SetParameters;

    CenterOfMassForceTransverser_(real4* force):force(force){}

    inline __device__ reduceType zero(const ComputationalData& computational,
                                      const SetParameters&    setParam,
                                      const int& reductionStep,
                                      const int& setIndex,
                                      const reduceType reduced[nReductions][nSet]){
        return make_real4(0.0);
    }

    inline __device__ reduceType transform(const int& index ,
                                           const ComputationalData& computational,
                                           const SetParameters&    setParam,
                                           const int& reductionStep,
                                           const int& setIndex,
                                           const reduceType reduced[nReductions][nSet]){

        const real m = computational.mass[index];
        return make_real4(m*make_real3(computational.pos[index]),m);
    }

    static inline __device__ reduceType reduce(const reduceType& reduced1,
                                               const reduceType& reduced2){
        return reduced1+reduced2;
    }

    //After reduction

    inline __device__ propertyType compute(const ComputationalData& computational,
                                           const SetParameters&    setParam,
                                           const reduceType reduced[nReductions][nSet]){

        const real  totalMass = reduced[0][0].w;
        const real3 com = computational.box.apply_pbc(make_real3(reduced[0][0])/totalMass);

        return PotentialType::force(com,computational,setParam);
    }

    inline __device__ spreadType spread(const propertyType& quantity,
                                        const int& index,
                                        const ComputationalData& computational,
                                        const SetParameters&    setParam,
                                        const int& setIndex,
                                        const reduceType reduced[nReductions][nSet]){

        const real totalMass = reduced[0][setIndex].w;

        return (computational.mass[index]/totalMass)*quantity;
    }

    inline __device__ void set(const int& index,
                               const spreadType& quantity,
                               const ComputationalData& computational,
                               const SetParameters&    setParam,
                               const int& setIndex,
                               const reduceType reduced[nReductions][nSet]){

        force[index] += quantity;
    }
};

template<class SetType_>
class CenterOfMass_{

    public:

        using SetType = SetType_;

        ///////////////////////////

        //Number of sets in the potential

        static constexpr int nSet = 1;
        static constexpr int nReductions = 1;

        std::vector<std::string> getSetLabels(){
            std::vector<std::string> labels = {"idSet_i"};
            return labels;
        }

        ///////////////////////////

        //Computational data
        using ComputationalData = typename SetType::ComputationalData;

        //Potential parameters
        using StorageData       = typename SetType::StorageData;

        //Set parameters
        using SetParameters     = typename SetType::SetParameters;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = CenterOfMassEnergyTransverser_<SetType,nSet,nReductions>;
        using ForceTransverser  = CenterOfMassForceTransverser_<SetType,nSet,nReductions>;

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return SetType::getComputationalData(this->gd,this->pd,storage,comp,st);
        }

        template<typename T>
        SetParameters processSetParameters(std::map<std::string,T>& setParametersMap){
            return SetType::processSetParameters(this->gd,setParametersMap);
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        CenterOfMass_(std::shared_ptr<GlobalData>  gd,
                      std::shared_ptr<ParticleGroup> pg,
                      DataEntry& data):gd(gd),
                                       pg(pg),pd(getExtendedParticleData(pg)){

            storage = SetType::getStorageData(gd,pd,data);
        }

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

}}}}
