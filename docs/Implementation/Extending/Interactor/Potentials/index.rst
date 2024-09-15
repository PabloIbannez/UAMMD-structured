Potentials
==========

UAMMD-structured uses a standard for Interactor templates. 
The class ``MyPotential_``, which expects to receive a template 
of Interactor ``MyType_<class Potential>``, must contain at least
two key attributes: ``StorageData`` and ``ComputationalData.`` 
Additionally, it must have the methods ``static __host__ getStorageData()`` 
and ``getComputationalData()`` to access the two attributes. 

``StorageData`` is intended to manage variables that will remain
constant throughout the simulation. For example, it can store
the parameters of a potential or the memory location for the
energy of the particles. The function ``getStorageData()``
is executed once in the simulation. On the other hand,
``ComputationalData`` contains the information sent to the GPU
to compute a computable. The function ``getComputationalData()``
is executed at each step of the simulation. Due to its nature,
any pointer contained in ComputationalData must refer to memory
stored on the GPU. Additionally, it is recommended that
ComputationalData be lightweight in terms of memory usage.

Finally, it should include ``static inline __device__ dataType computable()`` methods,
where the computation of different computables (energy, forces, etc.) 
due to the specific interaction is specified.

.. code-block:: cpp

    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "ParticleData/ExtendedParticleData.cuh"
    #include "ParticleData/ParticleGroup.cuh"
    
    #include "Interactor/Type/Type.cuh" //This file doesn't really exists
    #include "Interactor/InteractorFactory.cuh"
    
    namespace uammd{
    namespace structured{
    namespace Potentials{
    namespace Type{

        struct MyPotential_{
    
            struct StorageData{
                // Data storaged at the beggining of the simulation
                // GPU never use this information
            };

            struct ComputationalData{
                // Data required every step
                // GPU receive this information
            };


            //Storage data reader
            static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                       std::shared_ptr<ParticleGroup> pg,
                                                       DataEntry& data){
                StorageData storage;
                return storage;
            }

            //Computational data getter
            static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                                   std::shared_ptr<ParticleGroup> pg,
                                                                   const StorageData&  storage,
                                                                   const Computables& comp,
                                                                   const cudaStream_t& st){
    
                ComputationalData computational;
                // computational.example = storage.example
                return computational;
            }
    
            static inline __device__ real energy(int index_i, int index_j,
                                                 const ComputationalData& computational){
                real e = 0;
                // e is the energy of particle i due to particle j
                // associated wit this interaction
                return e;
    
            }
    
            static inline __device__ real3 force(int index_i, int index_j,
                                                 const ComputationalData& computational){

                real3 f = make_real3(0.0,0.0,0.0);
                // f is the force of particle i due to particle j
                // associated wit this interaction
                return f;
            }

            //static inline __device__ otherComputables() can be added here
    
        };
    
    }}}}
    
    REGISTER_TYPE_INTERACTOR( //replace TYPE for the actual Type
        Type,MyPotential,
        uammd::structured::Interactor::Type<MyPotential>
    )

It is highly recommended to review the specific template that
fits what you want to implement, as there are certain details
that differentiate one type of Interactor from another.
The templates presented here represent the minimum requirements
that a Potential must have to compile, but there is considerable
freedom to add intermediate classes or expand existing ones to meet the user's needs.

-------------------------------------------------

Available templates are:

.. toctree::
   :maxdepth: 2

   Bonds/index
   External
   Surface
   NonBonded

