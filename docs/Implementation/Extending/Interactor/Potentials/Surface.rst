Surface
========

Surface follows the minimal template example for potential. In this case, the
'data' and 'parameters' inputs do not have any mandatory information, so they
are entirely available to the user. Surface Interactor iterates over the
particles one by one, so the computables calculate properties for each particle.

.. code-block:: cpp

    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "ParticleData/ExtendedParticleData.cuh"
    #include "ParticleData/ParticleGroup.cuh"

    #include "Interactor/Single/SingleInteractor.cuh"
    #include "Interactor/Single/Surface/Surface.cuh"
    #include "Interactor/InteractorFactory.cuh"

    namespace uammd{
    namespace structured{
    namespace Potentials{
    namespace Surface{

        struct MyPotential_{

            // Computacional data
            struct ComputationalData{
                // Define data members specific to the potential
            };

            struct StorageData{
                // Define storage parameters specific to the potential
            };

            static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData> gd,
                                                                   std::shared_ptr<ParticleGroup> pg,
                                                                   const StorageData& storage,
                                                                   const Computables& comp,
                                                                   const cudaStream_t& st){

                ComputationalData computational;

                // Initialize computational data based on storage and global data

                return computational;
            }

            static __host__ StorageData getStorageData(std::shared_ptr<GlobalData> gd,
                                                       std::shared_ptr<ParticleGroup> pg,
                                                       DataEntry& data){

                StorageData storage;

                // Retrieve parameters from DataEntry and initialize StorageData

                return storage;
            }

            static inline __device__ real energy(int index_i, const ComputationalData& computational){
                // Implement energy computation specific to the potential
                return real(0.0);
            }

            static inline __device__ real3 force(int index_i, const ComputationalData& computational){
                // Implement force computation specific to the potential
                return make_real3(0.0);
            }
        };

        using MyPotential = Surface_<MyPotential_>;

    }}}}

    REGISTER_SINGLE_INTERACTOR(
        Surface,MyPotential,
        uammd::structured::Interactor::SingleInteractor<uammd::structured::Potentials::Surface::MyPotential>
    )

To register your own Surface potential create the file
``src/Interactor/Single/Surface/myPotential.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "Interactor":
        "Single":[
                 ["..."],
                 ["Surface","myPotential","myPotential.cu"]
                 ]
   }

