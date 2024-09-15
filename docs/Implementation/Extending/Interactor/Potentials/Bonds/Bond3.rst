Bond3
^^^^^

You can copy Bond2 example, substitute Bond2 with Bond3 in every part of the code, 
add ``int index_k`` to the computables functions arguments and create your 
generic Bond3 Potential. But most of the times, bond interactions between three particles 
only deppend on the angle between them. So UAMMD-structured offer the option
of create ``AngularBond3`` potentials.

The template is similar, but the only arguments of ``computable()`` functions are now
the angle between the three bonded particles and ComputationalData. Also, force() 
and hessian() are not defined. Instead you must define energyDerivate() and 
energySecondDerivate().

.. code-block:: cpp

    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "ParticleData/ExtendedParticleData.cuh"
    #include "ParticleData/ParticleGroup.cuh"

    #include "Interactor/Bonds/BondsInteractor.cuh"
    #include "Interactor/Bonds/Bond3/Bond3.cuh"
    #include "Interactor/InteractorFactory.cuh"

    namespace uammd{
    namespace structured{
    namespace Potentials{
    namespace Bond3{

        struct MyAngularPotential_{

            //Potential parameters

            struct ComputationalData{
                // Add necessary fields foor computational data
            };

            struct StorageData{
                // Add fields to store potential-related data
            };

            struct BondParameters{
                real3 parameter1; //Example, replace with actual parameters
            };

            //Computational data getter

            static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                                   std::shared_ptr<ParticleGroup> pg,
                                                                   const StorageData&  storage,
                                                                   const Computables& computables,
                                                                   const cudaStream_t&est){

                ComputationalData computational;

                return computational;
            }

            //Storage data reader

            static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                       std::shared_ptr<ParticleGroup> pg,
                                                       DataEntry& data){
                StorageData storage;
                return storage;
            }

            //Bond parameters reader

            template<typename T>
            static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                                 std::map<std::string,T>& bondParametersMap){

                BondParameters param;

                param.parameter1 = bondParametersMap.at("parameter1");

                return param;
            }

            //Energy and force definition

            static inline __device__ real energy(const real& ang,
                                                 const ComputationalData &computational,
                                                 const BondParameters &bondParam){

                const real energyValue = real(0.0);
                return energyValue;
            }

            static inline __device__ real energyDerivate(const real& ang,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){

                const real energyDerivateValue = real(0.0);
                return energyDerivateValue;
            }

            static inline __device__ real energySecondDerivate(const real& ang,
                                                               const ComputationalData &computational,
                                                               const BondParameters &bondParam){

                const real energySecondDerivateValue = real(0.0);
                return energySecondDerivateValue;
            }
        };

        // Alias the struct with AngularBond3 !!!
        using MyAngularPotential = AngularBond3_<MyAngularPotential_>;

    }}}}

    REGISTER_BOND_INTERACTOR(
            Bond3,MyAngularPotential,
            uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond3::MyAngularPotential>
            )


