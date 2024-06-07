#ifndef __MORSE_BOND2__
#define __MORSE_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct Morse_{

        struct ComputationalData{
            real4* pos;
            Box    box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{

            real r0;
            real E;
            real D;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box = gd->getEnsemble()->getBox();

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;
            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.r0 = bondParametersMap.at("r0");
            param.E  = bondParametersMap.at("E");
            param.D  = bondParametersMap.at("D");

            return param;
        }

        //Energy and force definition


        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real r0 = bondParam.r0;
            const real E  = bondParam.E;
            const real D  = bondParam.D;

            const real r2 = dot(rij, rij);

            real e = BasicPotentials::Morse::energy(rij,r2,E,r0,D);

            return e;
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real r0 = bondParam.r0;
            const real E  = bondParam.E;
            const real D  = bondParam.D;

            const real r2 = dot(rij, rij);

            real3 f = BasicPotentials::Morse::force(rij,r2,E,r0,D);

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;
        }


    };

    struct MorseCommon_D_{

        struct ComputationalData: public Morse_::ComputationalData{
            real D;
        };

        //Potential parameters

        struct StorageData: public Morse_::StorageData{
            real D;
        };

        struct BondParameters{

            real r0;
            real E;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;
            static_cast<Morse_::ComputationalData&>(computational)
            = Morse_::getComputationalData(gd,pg,storage,computables,st);

            computational.D = storage.D;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.D  = data.getParameter<real>("D");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.r0 = bondParametersMap.at("r0");
            param.E  = bondParametersMap.at("E");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Morse_::BondParameters bP;
            bP.D  = computational.D;
            bP.r0 = bondParam.r0;
            bP.E  = bondParam.E;

            return Morse_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Morse_::BondParameters bP;
            bP.D  = computational.D;
            bP.r0 = bondParam.r0;
            bP.E  = bondParam.E;

            return Morse_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

    };

    struct MorseCommon_r0_E_D_{

        struct ComputationalData: public Morse_::ComputationalData{
            real r0;
            real E;
            real D;
        };

        //Potential parameters

        struct StorageData: public Morse_::StorageData{
            real r0;
            real E;
            real D;
        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){
            ComputationalData computational;
            static_cast<Morse_::ComputationalData&>(computational)
            = Morse_::getComputationalData(gd,pg,storage,computables,st);

            computational.r0 = storage.r0;
            computational.E  = storage.E;
            computational.D  = storage.D;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                  std::shared_ptr<ParticleGroup> pg,
                                                  DataEntry& data){

            StorageData storage;

            storage.r0 = data.getParameter<real>("r0");
            storage.E  = data.getParameter<real>("E");
            storage.D  = data.getParameter<real>("D");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;
            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Morse_::BondParameters bP;
            bP.r0 = computational.r0;
            bP.E  = computational.E;
            bP.D  = computational.D;

            return Morse_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Morse_::BondParameters bP;
            bP.r0 = computational.r0;
            bP.E  = computational.E;
            bP.D  = computational.D;

            return Morse_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

    };

    struct MorseWCACommon_eps0_{

        struct ComputationalData{
            real4* pos;
            Box    box;
            real   eps0;
        };

        //Potential parameters

        struct StorageData{
            real eps0;
        };

        struct BondParameters{

            real r0;
            real E;
            real D;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box = gd->getEnsemble()->getBox();

            computational.eps0 = storage.eps0;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.eps0 = data.getParameter<real>("eps0");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.r0 = bondParametersMap.at("r0");
            param.E  = bondParametersMap.at("E");
            param.D  = bondParametersMap.at("D");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real eps0 = computational.eps0;

            const real r0 = bondParam.r0;
            const real E  = bondParam.E;
            const real D  = bondParam.D;

            const real r2 = dot(rij, rij);

            real e = BasicPotentials::Morse::energy(rij,r2,E,r0,D) +
                     BasicPotentials::WCA::Type2::energy(rij,r2,eps0,r0);

            return e;
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real eps0 = computational.eps0;

            const real r0 = bondParam.r0;
            const real E  = bondParam.E;
            const real D  = bondParam.D;

            const real r2 = dot(rij, rij);

            real3 f = BasicPotentials::Morse::force(rij,r2,E,r0,D) +
                      BasicPotentials::WCA::Type2::force(rij,r2,eps0,r0);

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;
        }

    };


    using Morse               = Bond2_<Morse_>;
    using MorseCommon_D       = Bond2_<MorseCommon_D_>;
    using MorseCommon_r0_E_D  = Bond2_<MorseCommon_r0_E_D_>;

    using MorseWCACommon_eps0 = Bond2_<MorseWCACommon_eps0_>;

}}}}

#endif
