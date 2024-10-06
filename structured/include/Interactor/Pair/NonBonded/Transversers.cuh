#pragma once

#include "uammd.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    //Simple transversers

    template <class NonBondedType_>
    struct EnergyTransverser_{

        real*  energy;

        using NonBondedType = NonBondedType_;
        using resultType    = real;

        EnergyTransverser_(std::shared_ptr<ParticleData> pd){
            this->energy = pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
        }

        inline __device__ resultType zero(){return real(0.0);}

        inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

        inline __device__ resultType compute(const int& index_i,
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){
            return NonBondedType::energy(index_i, index_j, computational)/real(2.0); //We divide by 2 because we are counting each interaction twice
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            energy[index_i] += quantity;
        }
    };

    template <class NonBondedType_>
    struct ForceTransverser_{

        real4*  force;

        using NonBondedType  = NonBondedType_;
        using resultType = real4;

        ForceTransverser_(std::shared_ptr<ParticleData> pd){
            this->force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
        }

        inline __device__ resultType zero(){return make_real4(0.0);}

        inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

        inline __device__ resultType compute(const int& index_i,
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){
            return make_real4(NonBondedType::force(index_i, index_j, computational),0.0);
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            force[index_i] += quantity;
        }
    };

    template <class NonBondedType_>
    struct LambdaTransverser_{

        real* lambdaDerivative;

        using NonBondedType  = NonBondedType_;
        using resultType = real;

        LambdaTransverser_(std::shared_ptr<ParticleData> pd){
            this->lambdaDerivative = pd->getLambdaDerivative(access::location::gpu, access::mode::readwrite).raw();
        }

        inline __device__ resultType zero(){return real(0.0);}

        inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

        inline __device__ resultType compute(const int& index_i,
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){
            return NonBondedType::lambdaDerivative(index_i, index_j, computational)/real(2.0); //We divide by 2 because we are counting each interaction twice
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            lambdaDerivative[index_i] += quantity;
        }
    };

    template <class NonBondedType_>
    struct StressTransverser_{

        tensor3*  stress;

        using NonBondedType  = NonBondedType_;
        using resultType = tensor3;

        StressTransverser_(std::shared_ptr<ParticleData> pd){
            this->stress = pd->getStress(access::location::gpu, access::mode::readwrite).raw();
        }

        inline __device__ resultType zero(){return tensor3(0.0);}

        inline __device__ void accumulate(resultType& total,const resultType current){
            total+=current;
        }

        inline __device__ resultType compute(const int& index_i,
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){
            real3 fij = NonBondedType::force(index_i, index_j, computational);
            real3 rij = computational.box.apply_pbc(make_real3(computational.pos[index_j]) - make_real3(computational.pos[index_i]));

            return computeStress(rij,fij);
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            stress[index_i] += quantity;
        }
    };

    template <class NonBondedType_>
    struct PairwiseForceTransverser_{

        real4*    pairwiseForce;
        const int*  id;
        const int*  selectedId;
        const int*  id2index;

        using NonBondedType   = NonBondedType_;
        using resultType = real4;

        PairwiseForceTransverser_(std::shared_ptr<ParticleData> pd){
            this->pairwiseForce = pd->getPairwiseForce(access::location::gpu, access::mode::readwrite).raw();
            this->id            = pd->getId(access::location::gpu, access::mode::read).raw();
            this->selectedId    = pd->getSelectedId(access::location::gpu, access::mode::read).raw();
            this->id2index      = pd->getIdOrderedIndices(access::location::gpu);
        }

        inline __device__ resultType zero(){return real4();}

        inline __device__ void accumulate(resultType& total,const resultType current){
            total+=current;
        }

        inline __device__ resultType compute(const int index_i,
                                             const int index_j,
                                             const typename NonBondedType::ComputationalData& computational){

            const int id_j  = id2index[index_j];
            const int selId = selectedId[index_i];

            real4 f = real4();
            if (id_j == selId){
                f = make_real4(NonBondedType::force(index_i,index_j,computational),0.0);
                //Force that particle j (selId here) makes over particle i
            }
            return f;
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            //I want to store the force that all the particles make over selId, but in compute I compute
            //the force that selId makes over the other particles.
            //That is the reason of the -= instead of +=
            pairwiseForce[index_i] -= quantity;
        }
    };

    template <class NonBondedType_>
    struct HessianTransverser_{

        tensor3*    hessian;
        const int*  id;
        const int*  selectedId;
        const int*  id2index;

        using NonBondedType   = NonBondedType_;
        using resultType      = tensor3;

        HessianTransverser_(std::shared_ptr<ParticleData> pd){
                this->hessian    = pd->getHessian(access::location::gpu, access::mode::readwrite).raw();
                this->id         = pd->getId(access::location::gpu, access::mode::read).raw();
                this->selectedId = pd->getSelectedId(access::location::gpu, access::mode::read).raw();
                this->id2index   = pd->getIdOrderedIndices(access::location::gpu);
        }

        inline __device__ resultType zero(){return tensor3(0.0);}

        inline __device__ void accumulate(resultType& total,const resultType current){
            total+=current;
        }

        inline __device__ resultType compute(const int& index_i,
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){

            // 11 12 13
            // 21 22 23
            // 31 32 33

            //We compute a column of the hessian
            //The element hessian[i] will store the block (i,selectedId)
            //We first derevite i and then selectedId

            const int id_i = id[index_i];
            const int id_j = id[index_j];

            const int selId     = selectedId[index_i];

            tensor3 H = tensor3(0.0);

            if (selId == id_i){ //Contribution to the self term
                                //self term has a minus sign
                H = -NonBondedType::hessian(id_i,id_j,computational);
            } else if (selId == id_j){
                H =  NonBondedType::hessian(id_i,id_j,computational);
            }

            return H;
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            hessian[index_i] += quantity;
        }
    };

    template <class NonBondedType_>
    struct MagneticFieldTransverser_{

        real4* magneticField;

        using NonBondedType = NonBondedType_;
        using resultType    = real4;

        MagneticFieldTransverser_(std::shared_ptr<ParticleData> pd){
            this->magneticField = pd->getMagneticField(access::location::gpu, access::mode::readwrite).raw();
        }

        inline __device__ resultType zero(){return make_real4(0.0);}

        inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

        inline __device__ resultType compute(const int& index_i,
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){
            return NonBondedType::magneticField(index_i, index_j, computational);
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            magneticField[index_i] += quantity;
        }
    };

    //Combined transversers

    template <class NonBondedType_>
    struct ForceTorqueTransverser_{

        real4*  force;
        real4*  torque;

        using NonBondedType = NonBondedType_;
        using resultType    = ForceTorque;

        ForceTorqueTransverser_(std::shared_ptr<ParticleData> pd){
            this->force  = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            this->torque = pd->getTorque(access::location::gpu, access::mode::readwrite).raw();
        }

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
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){
            return NonBondedType::forceTorque(index_i, index_j, computational);
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            force[index_i]  += quantity.force;
            torque[index_i] += quantity.torque;
        }
    };


    template <class NonBondedType_>
    struct ForceTorqueMagneticFieldTransverser_{

        real4*  force;
        real4*  torque;
        real4*  magneticField;

        using NonBondedType = NonBondedType_;
        using resultType    = ForceTorqueMagneticField;

        ForceTorqueMagneticFieldTransverser_(std::shared_ptr<ParticleData> pd){
            this->force         = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            this->torque        = pd->getTorque(access::location::gpu, access::mode::readwrite).raw();
            this->magneticField = pd->getMagneticField(access::location::gpu, access::mode::readwrite).raw();
        }

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
                                             const int& index_j,
                                             const typename NonBondedType::ComputationalData& computational){
            return NonBondedType::forceTorqueMagneticField(index_i, index_j, computational);
        }

        inline __device__ void set(const int& index_i,resultType& quantity){
            force[index_i]         += quantity.force;
            torque[index_i]        += quantity.torque;
            magneticField[index_i] += quantity.magneticField;
        }
    };

    #define DEFINE_TRANSVERSER_HAS_MEMBER(member_name)                     \
    template <typename T, typename = void>                                 \
    struct has_##member_name : std::false_type {};                         \
                                                                           \
    template <typename T>                                                  \
    struct has_##member_name<T, std::void_t<decltype(&T::member_name)>>    \
        : std::true_type {};

    DEFINE_TRANSVERSER_HAS_MEMBER(energy)
    DEFINE_TRANSVERSER_HAS_MEMBER(force)
    DEFINE_TRANSVERSER_HAS_MEMBER(lambdaDerivative)
    DEFINE_TRANSVERSER_HAS_MEMBER(stress)
    DEFINE_TRANSVERSER_HAS_MEMBER(hessian)
    DEFINE_TRANSVERSER_HAS_MEMBER(magneticField)
    DEFINE_TRANSVERSER_HAS_MEMBER(forceTorque)
    DEFINE_TRANSVERSER_HAS_MEMBER(forceTorqueMagneticField)

    #undef DEFINE_TRANSVERSER_HAS_MEMBER

}}}}
