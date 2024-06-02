#ifndef __GAUSSIAN_BOND2__
#define __GAUSSIAN_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct Gaussian_{

        struct ComputationalData{
            real4* pos;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{

            real E;
            real r0;
            real D;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos     = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box     = gd->getEnsemble()->getBox();

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

            param.E  = bondParametersMap.at("E");
            param.r0 = bondParametersMap.at("r0");
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

            const real E   = bondParam.E;
            const real r0  = bondParam.r0;
            const real D   = bondParam.D;

            const real r2 = dot(rij, rij);

            const real e = BasicPotentials::GaussianWell::energy(rij,r2,E,r0,D);

            return e;
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real E   = bondParam.E;
            const real r0  = bondParam.r0;
            const real D   = bondParam.D;

            const real r2 = dot(rij, rij);

            real3 f = BasicPotentials::GaussianWell::force(rij,r2,E,r0,D);

            if        (currentParticleIndex == index_i){

            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;
        }


    };

    struct GaussianCommon_E_r0_D_{

        struct ComputationalData : public Gaussian_::ComputationalData{
            real E;
            real r0;
            real D;
        };

        //Potential parameters

        struct StorageData: public Gaussian_::StorageData{
            real E;
            real r0;
            real D;
        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables){
            ComputationalData computational;

            static_cast<Gaussian_::ComputationalData&>(computational) = Gaussian_::getComputationalData(gd,pg,storage,computables);

            computational.E  = storage.E;
            computational.r0 = storage.r0;
            computational.D  = storage.D;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.E  = data.getParameter<real>("E");
            storage.r0 = data.getParameter<real>("r0");
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

            Gaussian_::BondParameters bP;
            bP.E  = computational.E;
            bP.r0 = computational.r0;
            bP.D  = computational.D;

            return Gaussian_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Gaussian_::BondParameters bP;
            bP.E  = computational.E;
            bP.r0 = computational.r0;
            bP.D  = computational.D;

            return Gaussian_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

    };

    using Gaussian              = Bond2_<Gaussian_>;
    using GaussianCommon_E_r0_D = Bond2_<GaussianCommon_E_r0_D_>;

}}}}

#endif
