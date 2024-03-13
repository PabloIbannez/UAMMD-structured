#ifndef __SET1_TORQUE_POTENTIAL__
#define __SET1_TORQUE_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Set1{

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
    using propertyType = real4;
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
            const real totalMass = reduced[0][setIndex].w;
            const real3 com = make_real3(reduced[0][setIndex])/totalMass;

            const real3 ri = make_real3(computational.pos[index])-com;

            const real3 T  = setParam.torque;
            const real3 ui = ri - T*dot(ri,T)/dot(T,T);

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

        return make_real4(0.0); //Not used
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

        //

        const real mi  = computational.mass[index];
        const real3 T  = setParam.torque;

        const real3 ri = make_real3(computational.pos[index])-com;
        const real3 ui = ri - T*dot(ri,T)/dot(T,T);

        const real3 f = mi*prefactor*cross(T,ui);

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

class Torque_ {

    public:

        ///////////////////////////

        //Number of sets in the potential

        static constexpr int nSet = 1;
        static constexpr int nReductions = 2;

        std::vector<std::string> getSetLabels(){
            std::vector<std::string> labels = {"idSet_i"};
            return labels;
        }

        ///////////////////////////

        //Computational data
        struct ComputationalData {
            real* mass;
            real4* pos;
        };

        //Set parameters
        struct SetParameters {
            real3 torque;
        };

        ///////////////////////////

        //Transversers
        using EnergyTransverser = TorqueEnergyTransverser_<ComputationalData,SetParameters,nSet,nReductions>;
        using ForceTransverser  = TorqueForceTransverser_<ComputationalData,SetParameters,nSet,nReductions>;

        ///////////////////////////

        ComputationalData getComputationalData(){
            ComputationalData computational;

            computational.pos  = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.mass = pd->getMass(access::location::gpu, access::mode::read).raw();

            return computational;
        }

        //Not storage needed

        template<typename T>
        SetParameters processSetParameters(std::map<std::string,T>& setParametersMap){

            SetParameters setParam;
            setParam.torque = setParametersMap["torque"];
            return setParam;
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

    public:

        Torque_(std::shared_ptr<GlobalData>    gd,
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

using ConstantTorqueOverCenterOfMass = Torque_;

}}}}

#endif
