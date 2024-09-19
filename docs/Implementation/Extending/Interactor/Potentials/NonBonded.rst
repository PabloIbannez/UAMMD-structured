NonBonded
=========

NonBonded is a type of Interactor designed to implement short-range interactions. 
It uses conditional neighbor lists, allowing you to define different interactions 
between different types of particles. The best part is that it is not different 
from implementing any other potential. You only need to define StorageData, 
ComputationalData, and the various computables(). The predefined NonBonded 
class manages the neighbor lists and efficiently iterates only over neighboring 
particles.

Remember to add the parameters cutOffFactor and condition to your parameter 
list in the input.json. Additionally, you should include which type of 
neighbourList you will use in the forceField section. For a more detailed 
explanation about how to set the input.json, refer to `ShortRange <../../../../Interactor/NonBonded/ShortRange/index.html>`_ .

.. code-block:: cpp

    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "ParticleData/ExtendedParticleData.cuh"
    #include "ParticleData/ParticleGroup.cuh"

    #include "Interactor/Pair/PairInteractor.cuh"
    #include "Interactor/Pair/NonBonded/NonBonded.cuh"
    #include "Interactor/InteractorFactory.cuh"

    #include "Interactor/BasicPotentials/myPotential.cuh" // This file don't really exists.

    namespace uammd{
    namespace structured{
    namespace Potentials{
    namespace NonBonded{

        struct myPotential_{

            //Computational data
            struct ComputationalData{

                real* pos;
                real cutOffFactor;
            };

            //Potential parameters
            struct StorageData{

                real cutOffFactor;
            };

            static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                                   std::shared_ptr<ParticleGroup> pg,
                                                                   const StorageData&  storage,
                                                                   const Computables& comp,
                                                                   const cudaStream_t& st){

                ComputationalData computational;

                std::shared_ptr<ParticleData> pd = pg->getParticleData();

                computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();

                computational.cutOffFactor = storage.cutOffFactor;

                return computational;
            }

            //Storage data reader

            static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                       std::shared_ptr<ParticleGroup> pg,
                                                       DataEntry& data){

                StorageData storage;

                storage.cutOffFactor = data.getParameter<real>("cutOffFactor");

                return storage;

            }


            static inline __device__ real energy(const int index_i,const int index_j,
                                                 const ComputationalData& computational){

                const real4 posi = computational.pos[index_i];
                const real4 posj = computational.pos[index_j];

                const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                real e = real(0.0);

                }

                return e;
            }

          static inline __device__ real3 force(const int index_i,const int index_j,
                                                 const ComputationalData& computational){

                const real4 posi = computational.pos[index_i];
                const real4 posj = computational.pos[index_j];

                const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                real3 f = make_real3(real(0.0));

                }

                return f;
            }

          static inline __device__ tensor3 hessian(const int index_i,const int index_j,
                               const ComputationalData& computational){

                const real4 posi = computational.pos[index_i];
                const real4 posj = computational.pos[index_j];

                const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                tensor3 H = tensor3(real(0.0));

                }

                return H;
            }
        };

        using myPotential = NonBondedHessian_<myPotential_>;

    }}}}

    REGISTER_NONBONDED_INTERACTOR(
        NonBonded,myPotential,
        uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::myPotential>
    )

To register your own NonBonded potential create the file
``src/Interactor/Pair/NonBonded/myPotential.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "Interactor":
        "Pair":[
                 ["..."],
                 ["NonBonded","myPotential","myPotential.cu"]
                 ]
   }


