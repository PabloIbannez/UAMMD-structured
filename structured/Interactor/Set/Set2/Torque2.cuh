#ifndef __SET2_TORQUE_POTENTIAL__
#define __SET2_TORQUE_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Set2{

template <class ComputationalData,class SetParameters,
          int nSet, int nReductions>
struct TorqueEnergyTransverser_{

    real*  energy;

    using reduceType   = real4;
    using propertyType = real;
    using spreadType   = real;

    TorqueEnergyTransverser_(real*  energy):energy(energy){}

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

        return make_real4(0.0);
    }

    static inline __device__ reduceType reduce(const reduceType& reduced1,
                                               const reduceType& reduced2){
        return reduced1+reduced2;
    }

    //After reduction

    inline __device__ propertyType compute(const ComputationalData& computational,
                                           const SetParameters&    setParam,
                                           const reduceType reduced[nReductions][nSet]){

        return real(0.0);
    }

    inline __device__ spreadType spread(const propertyType& quantity,
                                        const int& index,
                                        const ComputationalData& computational,
                                        const SetParameters&    setParam,
                                        const int& setIndex,
                                        const reduceType reduced[nReductions][nSet]){

        return quantity;
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

template <class ComputationalData,class SetParameters,
          int nSet, int nReductions>
struct TorqueForceTransverser_{

    real4* force;

    using reduceType   = real4;
    using propertyType = real3;
    using spreadType   = real4;

    TorqueForceTransverser_(real4*  force):force(force){}

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

        reduceType acc = make_real4(0.0);

        const real m = computational.mass[index];

        if (reductionStep == 0){
            acc = make_real4(m*make_real3(computational.pos[index]),m);
        }

        if (reductionStep == 1){

            //Using previus step
            const real totalMass0 = reduced[0][0].w;
            const real3 com0      = make_real3(reduced[0][0])/totalMass0;

            const real totalMass1 = reduced[0][1].w;
            const real3 com1      = make_real3(reduced[0][1])/totalMass1;

            real3 torqueVector = setParam.torque*normalize(computational.box.apply_pbc(com1-com0));
            if(setIndex == 1){
                torqueVector = -torqueVector;
            }

            real3 ri;
            if(setIndex == 0){
                ri = make_real3(computational.pos[index])-com0;
            } else {
                ri = make_real3(computational.pos[index])-com1;
            }

            const real3 ui = ri-torqueVector*dot(ri,torqueVector)/dot(torqueVector,torqueVector);

            acc.x = m*dot(ui,ui); // acc = (1/prefactor_i,0,0,0)
        }

        return acc;
    }

    static inline __device__ reduceType reduce(const reduceType& reduced1,
                                               const reduceType& reduced2){
        return reduced1+reduced2;
    }

    //After reduction

    inline __device__ propertyType compute(const ComputationalData& computational,
                                           const SetParameters&    setParam,
                                           const reduceType reduced[nReductions][nSet]){

        const real totalMass0 = reduced[0][0].w;
        const real3 com0      = make_real3(reduced[0][0])/totalMass0;

        const real totalMass1 = reduced[0][1].w;
        const real3 com1      = make_real3(reduced[0][1])/totalMass1;

        real3 torqueVector = setParam.torque*normalize(computational.box.apply_pbc(com1-com0));

        return torqueVector;
    }

    inline __device__ spreadType spread(const propertyType& quantity,
                                        const int& index,
                                        const ComputationalData& computational,
                                        const SetParameters&    setParam,
                                        const int& setIndex,
                                        const reduceType reduced[nReductions][nSet]){

        spreadType result;

        const real totalMass = reduced[0][setIndex].w;
        const real3 com      = make_real3(reduced[0][setIndex])/totalMass;

        const real prefactor = real(1.0)/reduced[1][setIndex].x;

        real3 torqueVector   = quantity;
        if(setIndex == 1){
            torqueVector = -torqueVector;
        }

        const real m = computational.mass[index];

        const real3 ri = make_real3(computational.pos[index])-com;
        const real3 ui = ri-torqueVector*dot(ri,torqueVector)/dot(torqueVector,torqueVector);

        const real3 f =  m*prefactor*cross(torqueVector,ui);

        return make_real4(f,0.0);
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

class Torque_{

    public:

        ///////////////////////////

        //Number of sets in the potential

        static constexpr int nSet = 2;
        static constexpr int nReductions = 2;

        std::vector<std::string> getSetLabels(){
            std::vector<std::string> labels = {"idSet_i","idSet_j"};
            return labels;
        }

        ///////////////////////////

        //Computational data
        struct ComputationalData {
            real* mass;
            real4* pos;
            Box box;
        };

        //Not storage needed

        //Set parameters
        struct SetParameters {
            real torque;
        };

        ///////////////////////////

        //Transversers
        using EnergyTransverser = TorqueEnergyTransverser_<ComputationalData,SetParameters,nSet,nReductions>;
        using ForceTransverser  = TorqueForceTransverser_<ComputationalData,SetParameters,nSet,nReductions>;

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp){
            ComputationalData computational;

            computational.pos  = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.mass = pd->getMass(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            return computational;
        }

        template<typename T>
        SetParameters processSetParameters(std::map<std::string,T>& setParametersMap){
            SetParameters setParam;
            setParam.torque = setParametersMap.at("torque");
            return setParam;
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

    public:

        Torque_(std::shared_ptr<GlobalData>  gd,
                std::shared_ptr<ParticleGroup> pg,
                DataEntry& data):gd(gd),
                                 pg(pg),pd(getExtendedParticleData(pg)){}


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

using ConstantTorqueBetweenCentersOfMass  = Torque_;

}}}}

#endif
