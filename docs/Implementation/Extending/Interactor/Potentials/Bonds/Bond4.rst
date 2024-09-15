Bond4
^^^^^

You can copy Bond2 example, substitute Bond2 with Bond4 in every part of the code,
add ``int index_k`` and ``int index_w`` to the computables functions arguments and create your
generic Bond4 Potential. But most of the times, bond interactions between four particles
only deppend on the dihedral angle between them, see `Bond4 <../../../../Interactor/Bonds/Bond4/index.html>`_. So UAMMD-structured offer the option
of create ``AngularBond4`` potentials.

The template is similar, but now ``computable()`` functions receive as arguments
:math:`sin(\phi)` and :math:`cos(\phi)` where :math:`\phi` is the dihedral angle
and ComputationalData. The decision to use the angle projections instead
of the angle value itself is because obtaining the projections requires
one less operation, and in many cases, the potentials depend on the projections.
You always can obtain the dihedral angle using :math:`\phi = atan2(sin(\phi),cos(\phi))`

Also, force() and hessian() are not defined. Instead you must define energyDerivate() and
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

            static inline __device__ real energy(const real& cos_dih, const real& sin_dih,
                                                 const ComputationalData &computational,
                                                 const BondParameters &bondParam){

                const real energyValue = real(0.0);
                return energyValue;
            }

            static inline __device__ real energyDerivate(const real& cos_dih, const real& sin_dih,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){

                const real energyDerivateValue = real(0.0);
                return energyDerivateValue;
            }

            static inline __device__ real energySecondDerivate(const real& cos_dih, const real& sin_dih,
                                                               const ComputationalData &computational,
                                                               const BondParameters &bondParam){

                const real energySecondDerivateValue = real(0.0);
                return energySecondDerivateValue;
            }
        };

        // Alias the struct with AngularBond3 !!!
        using MyAngularPotential = AngularBond4_<MyAngularPotential_>;

    }}}}

    REGISTER_BOND_INTERACTOR(
            Bond4,MyAngularPotential,
            uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond4::MyAngularPotential>
            )


