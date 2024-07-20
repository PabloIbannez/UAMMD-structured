#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    template <class NonBondedType_>
    struct EnergyTransverser_{

        real*  energy;

        using NonBondedType = NonBondedType_;
        using resultType    = real;

        EnergyTransverser_(real*  energy):energy(energy){}

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

        ForceTransverser_(real4*  force):force(force){}

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

        LambdaTransverser_(real* lambdaDerivative):lambdaDerivative(lambdaDerivative){}

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

        StressTransverser_(tensor3*  stress):stress(stress){}

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
    struct HessianTransverser_{

        tensor3*    hessian;
        const int*  id;
        const int*  selectedId;
        const int*  id2index;

        using NonBondedType   = NonBondedType_;
        using resultType      = tensor3;

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
    struct PairwiseForceTransverser_{

        real4*    pairwiseForce;
        const int*  id;
        const int*  selectedId;
        const int*  id2index;

        using NonBondedType   = NonBondedType_;
        using resultType = real4;

        PairwiseForceTransverser_(real4* pairwiseForce,
                const int* id,
                const int* selectedId,
                const int* id2index):pairwiseForce(pairwiseForce),
        id(id),
        selectedId(selectedId),
        id2index(id2index){}

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
    struct ForceTorqueTransverser_{

        real4*  force;
        real4*  torque;

        using NonBondedType = NonBondedType_;
        using resultType    = ForceTorque;

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
    struct MagneticFieldTransverser_{

        real4* magneticField;

        using NonBondedType = NonBondedType_;
        using resultType    = real4;

        MagneticFieldTransverser_(real4* magneticField):magneticField(magneticField){}

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

    template <class NonBondedType_>
    struct ForceTorqueMagneticFieldTransverser_{

        real4*  force;
        real4*  torque;
        real4*  magneticField;

        using NonBondedType = NonBondedType_;
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
                const int& index_j,
                const typename NonBondedType::ComputationalData& computational){
            ForceTorque forceTorque = NonBondedType::forceTorque(index_i, index_j, computational);
            real4 magneticField     = NonBondedType::magneticField(index_i, index_j, computational);

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

    ///////////////////////////////////////////////////////////////

    template<class NonBondedType_>
    class NonBondedBase_{

        public:

            using NonBondedType = NonBondedType_;

            //Computational data
            using ComputationalData = typename NonBondedType::ComputationalData;

            //Potential parameters
            using StorageData       = typename NonBondedType::StorageData;

            ///////////////////////////

            ComputationalData getComputationalData(const Computables& comp,
                    const cudaStream_t& st){
                return NonBondedType::getComputationalData(this->gd,this->pg,storage,comp,st);
            }

        protected:

            std::shared_ptr<GlobalData>    gd;
            std::shared_ptr<ParticleGroup> pg;

            std::shared_ptr<ExtendedParticleData> pd;

            StorageData storage;

        public:

            NonBondedBase_(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data):gd(gd),
            pg(pg),pd(getExtendedParticleData(pg)){

                storage = NonBondedType::getStorageData(gd,pg,data);

            }

            ///////////////////////////

            real getCutOff(){return storage.cutOff;}
    };

    template<class NonBondedType_>
    class NonBonded_ : public NonBondedBase_<NonBondedType_>{

        public:

            using NonBondedType = typename NonBondedBase_<NonBondedType_>::NonBondedType;

            ///////////////////////////

            //Transversers

            using EnergyTransverser = EnergyTransverser_<NonBondedType>;
            using ForceTransverser  = ForceTransverser_<NonBondedType>;
            using StressTransverser = StressTransverser_<NonBondedType>;
            using PairwiseForceTransverser = PairwiseForceTransverser_<NonBondedType>;

            ///////////////////////////

            NonBonded_(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data):NonBondedBase_<NonBondedType_>(gd,pg,data){}

            ///////////////////////////

            EnergyTransverser getEnergyTransverser(){

                real*  energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();

                return EnergyTransverser(energy);
            }

            ForceTransverser getForceTransverser(){

                real4*  force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                return ForceTransverser(force);
            }

            StressTransverser getStressTransverser(){

                tensor3*  stress = this->pd->getStress(access::location::gpu, access::mode::readwrite).raw();

                return StressTransverser(stress);
            }

            PairwiseForceTransverser getPairwiseForceTransverser(){

                real4*  pforce        = this->pd->getPairwiseForce(access::location::gpu, access::mode::readwrite).raw();
                const int* id2index   = this->pd->getIdOrderedIndices(access::location::gpu);
                const int* id         = this->pd->getId(access::location::gpu, access::mode::read).raw();
                const int* selectedId = this->pd->getSelectedId(access::location::gpu, access::mode::read).raw();


                return PairwiseForceTransverser(pforce,
                        id,
                        selectedId,id2index);
            }
    };

    template<class NonBondedType_>
    class NonBondedLambda_ : public NonBonded_<NonBondedType_>{

        public:

            using NonBondedType = typename NonBonded_<NonBondedType_>::NonBondedType;

            ///////////////////////////

            //Transversers

            using LambdaTransverser = LambdaTransverser_<NonBondedType>;

            ///////////////////////////

            NonBondedLambda_(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data):NonBonded_<NonBondedType_>(gd,pg,data){}

            ///////////////////////////

            LambdaTransverser getLambdaTransverser(){

                real*  lambdaDerivative = this->pd->getLambdaDerivative(access::location::gpu, access::mode::readwrite).raw();

                return LambdaTransverser(lambdaDerivative);
            }
    };

    template<class NonBondedType_>
    class NonBondedTorque_ : public NonBondedBase_<NonBondedType_>{

        public:

            using NonBondedType = typename NonBondedBase_<NonBondedType_>::NonBondedType;

            ///////////////////////////

            //Transversers

            using EnergyTransverser = EnergyTransverser_<NonBondedType>;
            using ForceTransverser  = ForceTorqueTransverser_<NonBondedType>;

            ///////////////////////////

            NonBondedTorque_(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data):NonBondedBase_<NonBondedType_>(gd,pg,data){}

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

    template<class NonBondedType_>
    class NonBondedForceTorqueMagneticField_ : public NonBondedTorque_<NonBondedType_>{

        public:

            using NonBondedType = typename NonBondedTorque_<NonBondedType_>::NonBondedType;

            ///////////////////////////

            //Transversers

            using MagneticFieldTransverser = MagneticFieldTransverser_<NonBondedType>;
            using ForceTorqueMagneticFieldTransverser = ForceTorqueMagneticFieldTransverser_<NonBondedType>;

            ///////////////////////////

            NonBondedForceTorqueMagneticField_(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data):NonBondedTorque_<NonBondedType_>(gd,pg,data){}

            ///////////////////////////

            MagneticFieldTransverser getMagneticFieldTransverser(){

                real4*  magneticField = this->pd->getMagneticField(access::location::gpu, access::mode::readwrite).raw();

                return MagneticFieldTransverser(magneticField);
            }

            ForceTorqueMagneticFieldTransverser getForceTorqueMagneticFieldTransverser(){

                real4*  magneticField = this->pd->getMagneticField(access::location::gpu, access::mode::readwrite).raw();
                real4*  force         = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
                real4*  torque        = this->pd->getTorque(access::location::gpu, access::mode::readwrite).raw();
                return ForceTorqueMagneticFieldTransverser(force, torque, magneticField);
            }
    };

    template<class NonBondedType_>
    class NonBondedHessian_ : public NonBonded_<NonBondedType_> {

        public:

            using NonBondedType = typename NonBonded_<NonBondedType_>::NonBondedType;

            ///////////////////////////

            //Transverser
            using HessianTransverser = HessianTransverser_<NonBondedType>;

        public:

            NonBondedHessian_(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data):NonBonded_<NonBondedType_>(gd,pg,data){}

            ///////////////////////////

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

}}}}
