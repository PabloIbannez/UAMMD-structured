#ifndef __CONSTANT_FORCE_BOND1__
#define __CONSTANT_FORCE_BOND1__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond1{

    struct ConstantForce_{

        struct ComputationalData{};

        struct StorageData{};

        struct BondParameters{
            real3 constantForce;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

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

            param.constantForce    = bondParametersMap.at("force");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters    &bondParam){

            const real e = real(0.0);
            return e;
        }


        static inline __device__ real3 force(int index_i,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters    &bondParam){

            const real3 f = make_real3(bondParam.constantForce);

            return f;
        }

    };

    //We rename the struct to make it easier to use
    using ConstantForce = Bond1_<ConstantForce_>;

}}}}

#endif
