#pragma once

#include "Interactor/Pair/NonBonded/Transversers.cuh"

#define DECLARE_TRANSVERSER(name, Name) \
    template<typename T, typename = void> \
    struct Name##Transverser { \
        using impl = void; \
    }; \
    \
    template<typename T> \
    struct Name##Transverser<T, std::enable_if_t<has_##name<T>::value>> { \
        using impl = Name##Transverser_<T>; \
    };

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    template<class NonBondedType_>
    class NonBonded_{

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

            ///////////////////////////

            // Conditionally define all transverser types

            struct Transversers_{
                // Simple transversers
                DECLARE_TRANSVERSER(energy, Energy)
                DECLARE_TRANSVERSER(force, Force)
                DECLARE_TRANSVERSER(lambdaDerivative, Lambda)
                DECLARE_TRANSVERSER(stress, Stress)
                DECLARE_TRANSVERSER(force, PairwiseForce)
                DECLARE_TRANSVERSER(hessian, Hessian)
                DECLARE_TRANSVERSER(magneticField, MagneticField)

                // Combined transversers
                DECLARE_TRANSVERSER(forceTorque, ForceTorque)
                DECLARE_TRANSVERSER(forceTorqueMagneticField, ForceTorqueMagneticField)
            };

            // Simple transversers
            using EnergyTransverser        = typename Transversers_::EnergyTransverser<NonBondedType>::impl;
            using ForceTransverser         = typename Transversers_::ForceTransverser<NonBondedType>::impl;
            using LambdaTransverser        = typename Transversers_::LambdaTransverser<NonBondedType>::impl;
            using StressTransverser        = typename Transversers_::StressTransverser<NonBondedType>::impl;
            using PairwiseForceTransverser = typename Transversers_::PairwiseForceTransverser<NonBondedType>::impl;
            using HessianTransverser       = typename Transversers_::HessianTransverser<NonBondedType>::impl;
            using MagneticFieldTransverser = typename Transversers_::MagneticFieldTransverser<NonBondedType>::impl;

            // Combined transversers
            using ForceTorqueTransverser              = typename Transversers_::ForceTorqueTransverser<NonBondedType>::impl;
            using ForceTorqueMagneticFieldTransverser = typename Transversers_::ForceTorqueMagneticFieldTransverser<NonBondedType>::impl;

            // Check if transversers are available
            template<typename T>
            struct is_transverser_available : std::bool_constant<!std::is_same_v<T, void>> {};

            // Simple transversers
            static constexpr bool isEnergyTransverserAvailable        = is_transverser_available<EnergyTransverser>::value;
            static constexpr bool isForceTransverserAvailable         = is_transverser_available<ForceTransverser>::value;
            static constexpr bool isLambdaTransverserAvailable        = is_transverser_available<LambdaTransverser>::value;
            static constexpr bool isStressTransverserAvailable        = is_transverser_available<StressTransverser>::value;
            static constexpr bool isPairwiseForceTransverserAvailable = is_transverser_available<PairwiseForceTransverser>::value;
            static constexpr bool isHessianTransverserAvailable       = is_transverser_available<HessianTransverser>::value;
            static constexpr bool isMagneticFieldTransverserAvailable = is_transverser_available<MagneticFieldTransverser>::value;

            // Combined transversers
            static constexpr bool isForceTorqueTransverserAvailable              = is_transverser_available<ForceTorqueTransverser>::value;
            static constexpr bool isForceTorqueMagneticFieldTransverserAvailable = is_transverser_available<ForceTorqueMagneticFieldTransverser>::value;


        protected:

            std::shared_ptr<GlobalData>    gd;
            std::shared_ptr<ParticleGroup> pg;

            std::shared_ptr<ExtendedParticleData> pd;

            StorageData storage;

        public:

            NonBonded_(std::shared_ptr<GlobalData>    gd,
                       std::shared_ptr<ParticleGroup> pg,
                       DataEntry& data):
                       gd(gd),pg(pg),pd(getExtendedParticleData(pg)){
                storage = NonBondedType::getStorageData(gd,pg,data);
            }

            ///////////////////////////

            real getCutOff(){return storage.cutOff;}

            ///////////////////////////

            // Getters for the transversers

            // Simple transversers
            EnergyTransverser getEnergyTransverser() {
                if constexpr (isEnergyTransverserAvailable) {
                    return EnergyTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Energy transverser not available.");
                    return EnergyTransverser();
                }
            }

            ForceTransverser getForceTransverser() {
                if constexpr (isForceTransverserAvailable) {
                    return ForceTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Force transverser not available.");
                    return ForceTransverser();
                }
            }

            LambdaTransverser getLambdaTransverser() {
                if constexpr (isLambdaTransverserAvailable) {
                    return LambdaTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Lambda transverser not available.");
                    return LambdaTransverser();
                }
            }

            StressTransverser getStressTransverser() {
                if constexpr (isStressTransverserAvailable) {
                    return StressTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Stress transverser not available.");
                    return StressTransverser();
                }
            }

            PairwiseForceTransverser getPairwiseForceTransverser() {
                if constexpr (isPairwiseForceTransverserAvailable) {
                    return PairwiseForceTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Pairwise force transverser not available.");
                    return PairwiseForceTransverser();
                }
            }

            HessianTransverser getHessianTransverser() {
                if constexpr (isHessianTransverserAvailable) {
                    return HessianTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Hessian transverser not available.");
                    return HessianTransverser();
                }
            }

            MagneticFieldTransverser getMagneticFieldTransverser() {
                if constexpr (isMagneticFieldTransverserAvailable) {
                    return MagneticFieldTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Magnetic field transverser not available.");
                    return MagneticFieldTransverser();
                }
            }

            // Combined transversers
            ForceTorqueTransverser getForceTorqueTransverser() {
                if constexpr (isForceTorqueTransverserAvailable) {
                    return ForceTorqueTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Force-torque transverser not available.");
                    return ForceTorqueTransverser();
                }
            }

            ForceTorqueMagneticFieldTransverser getForceTorqueMagneticFieldTransverser() {
                if constexpr (isForceTorqueMagneticFieldTransverserAvailable) {
                    return ForceTorqueMagneticFieldTransverser(this->pd);
                } else {
                    System::log<System::CRITICAL>("[NonBonded] Force-torque-magnetic field transverser not available.");
                    return ForceTorqueMagneticFieldTransverser();
                }
            }

    };

}}}}
