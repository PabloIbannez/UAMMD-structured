Bond1
^^^^^

Iterates over every particle for each Bond, so computable functions expect to receive
particle index (``currentParticleIndex``) and bond index (``index_i``).

.. code-block:: cpp

    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "ParticleData/ExtendedParticleData.cuh"
    #include "ParticleData/ParticleGroup.cuh"

    #include "Interactor/Bonds/BondsInteractor.cuh"
    #include "Interactor/Bonds/Bond1/Bond1.cuh"
    #include "Interactor/InteractorFactory.cuh"

    #include "Interactor/BasicPotentials/myPotential.cuh"  // Replace with the actual potential file
    // if your creating a new BasicPotential, remember to add it to Interactor/BasicPotentials to use
    // it in other places of the code.

    namespace uammd {
    namespace structured {
    namespace Potentials {
    namespace Bond1 {

        struct myPotential_ {

            struct ComputationalData {
                real4* pos; // In case you are using positions to calculate computables.
                // Add necessary fields for computational data
            };

            struct StorageData {
                // Add fields to store potential-related data
            };

            struct BondParameters {
                real3 parameter1;  // Example, replace with actual parameters
            };

            // Computational data getter
            static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData> gd,
                                                                   std::shared_ptr<ParticleGroup> pg,
                                                                   const StorageData& storage,
                                                                   const Computables& computables,
                                                                   const cudaStream_t& st) {

                ComputationalData computational;

                std::shared_ptr<ParticleData> pd = pg->getParticleData();
                computational.pos = pd->getPos(access::location::gpu,access::mode::read).raw(); //You probably need positions to calculate forces or energy, if not delete this line
                // Fill computational data
                return computational;
            }

            // Storage data reader
            static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                       std::shared_ptr<ParticleGroup> pg,
                                                       DataEntry& data) {

                StorageData storage;
                // Initialize storage data fields
                return storage;
            }

            // Bond parameters reader
            template<typename T>
            static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                                 std::map<std::string, T>& bondParametersMap) {

                BondParameters param;
                param.parameter1 = bondParametersMap.at("parameter1");  // Example parameter
                return param;
            }

            // Force calculation
            static inline __device__ real3 force(int index_i,
                                                 int currentParticleIndex,
                                                 const ComputationalData &computational,
                                                 const BondParameters &bondParam) {

                const real3 posi = make_real3(computational.pos[index_i]);
                // Replace this with actual force calculation logic
                real3 forceValue = make_real3(0.0f, 0.0f, 0.0f);
                return forceValue;
            }

            // Energy calculation
            static inline __device__ real energy(int index_i,
                                                 int currentParticleIndex,
                                                 const ComputationalData &computational,
                                                 const BondParameters &bondParam) {

                const real3 posi = make_real3(computational.pos[index_i]);
                // Replace this with actual energy calculation logic
                const real energyValue = real(0.0);
                return energyValue;
            }

            /* static inline __device__ dataType computable(int index_i,
                                                            int currentParticleIndex,
                                                            const ComputationalData &computational,
                                                            const BondParameters &bondParam) {

                dataType C = dataType();
                return C;
            }*/

        };

        // Alias the struct for ease of use
        using myPotential = Bond1_<myPotential_>;

    }}}}

    REGISTER_BOND_INTERACTOR(
        Bond1, myPotential,
        uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond1::myPotential>
    )

To register your own Bond1 potential create the file
``src/Interactor/Bonds/Bond1/myPotential.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "Interactor":
        "Bonds":[
                ["..."],
                ["Bond1","myPotential","myPotential.cu"]
                ]
   }

