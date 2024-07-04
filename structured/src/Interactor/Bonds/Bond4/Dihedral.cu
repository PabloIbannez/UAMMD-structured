#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond4/Bond4.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond4{

    struct Dihedral_{

        struct ComputationalData{
            real4* pos;
            Box    box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{
            int  n;
            real K;
            real phi0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.box = gd->getEnsemble()->getBox();
            computational.pos = pd->getPos(access::location::gpu,
                                           access::mode::read).raw();

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

            param.n    = bondParametersMap.at("n");
            param.K    = bondParametersMap.at("K");
            param.phi0 = bondParametersMap.at("phi0");

            return param;
        }

        //Angular potential functions

        static inline __device__ real energy(const real& cos_dih,const real& sin_dih,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            const int  n   = bondParam.n;
            const real K   = bondParam.K;
            const real pha = bondParam.phi0;

            real cosnt=real(1);
            real sinnt=real(0);

            real tmp;

            for(int krot=0;krot<n;krot++){
                tmp   = cosnt*cos_dih - sinnt*sin_dih;
                sinnt = sinnt*cos_dih + cosnt*sin_dih;
                cosnt = tmp;
            }

            const real cospha = cos(pha);
            const real sinpha = sin(pha);

            //K(1.0+cos(n*phi+phi_0))
            return K*(real(1.0)+ cospha*cosnt + sinnt*sinpha);
        }

        static inline __device__ real energyDerivate(const real& cos_dih,const real& sin_dih,
                                                     const ComputationalData &computational,
                                                     const BondParameters &bondParam){

            const int  n   = bondParam.n;
            const real K   = bondParam.K;
            const real pha = bondParam.phi0;

            real cosnt=real(1);
            real sinnt=real(0);

            //sinnt cosnt computation
            real tmp;

            for(int krot=0;krot<n;krot++){
                tmp   = cosnt*cos_dih - sinnt*sin_dih;
                sinnt = sinnt*cos_dih + cosnt*sin_dih;
                cosnt = tmp;
            }

            const real cospha = cos(pha);
            const real sinpha = sin(pha);

            return -K*n*(cospha*sinnt-cosnt*sinpha);
        }

      static inline __device__ real energySecondDerivate(const real& cos_dih,const real& sin_dih,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){

            const int  n   = bondParam.n;
            const real K   = bondParam.K;
            const real pha = bondParam.phi0;

            real cosnt=real(1);
            real sinnt=real(0);

            //sinnt cosnt computation
            real tmp;

            for(int krot=0;krot<n;krot++){
                tmp   = cosnt*cos_dih - sinnt*sin_dih;
                sinnt = sinnt*cos_dih + cosnt*sin_dih;
                cosnt = tmp;
            }

            const real cospha = cos(pha);
            const real sinpha = sin(pha);

            return -K*n*n*(cospha*cosnt+sinnt*sinpha);
        }

    };

    struct DihedralCommon_n_K_phi0_{

        struct ComputationalData: public Dihedral_::ComputationalData{
            int n;
            real K;
            real phi0;
        };

        //Potential parameters

        struct StorageData: public Dihedral_::StorageData{

            int  n;
            real K;
            real phi0;

        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){
            ComputationalData computational;
            static_cast<Dihedral_::ComputationalData&>(computational) =
            Dihedral_::getComputationalData(gd,pg,storage,computables,st);

            computational.n    = storage.n;
            computational.K    = storage.K;
            computational.phi0 = storage.phi0;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.n    = data.getParameter<int>("n");
            storage.K    = data.getParameter<real>("K");
            storage.phi0 = data.getParameter<real>("phi0");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;
            return param;
        }

        //Angular potential functions

        static inline __device__ real energy(const real& cos_dih,const real& sin_dih,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            Dihedral_::BondParameters bP;
            bP.n    = computational.n;
            bP.K    = computational.K;
            bP.phi0 = computational.phi0;

            return Dihedral_::energy(cos_dih,sin_dih,computational,bP);
        }

        static inline __device__ real energyDerivate(const real& cos_dih,const real& sin_dih,
                                                     const ComputationalData &computational,
                                                     const BondParameters &bondParam){

            Dihedral_::BondParameters bP;
            bP.n    = computational.n;
            bP.K    = computational.K;
            bP.phi0 = computational.phi0;

            return Dihedral_::energyDerivate(cos_dih,sin_dih,computational,bP);
        }

      static inline __device__ real energySecondDerivate(const real& cos_dih,const real& sin_dih,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){

            Dihedral_::BondParameters bP;
            bP.n    = computational.n;
            bP.K    = computational.K;
            bP.phi0 = computational.phi0;

            return Dihedral_::energySecondDerivate(cos_dih,sin_dih,computational,bP);
        }
    };

    struct Dihedral4_ {

        //Potential parameters

        using ComputationalData = Dihedral_::ComputationalData;
        using StorageData = Dihedral_::StorageData;

        struct BondParameters{

            real4 K;
            real4 phi0;

        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            return Dihedral_::getComputationalData(gd,pg,
                                                   storage,
                                                   computables,st);
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
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

            param.K    = bondParametersMap.at("K");
            param.phi0 = bondParametersMap.at("phi0");

            return param;
        }
        //Angular potential functions

        static inline __device__ real energy(const real& cos_dih,const real& sin_dih,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real K[4]  ={bondParam.K.x,
                         bondParam.K.y,
                         bondParam.K.z,
                         bondParam.K.w};

            real pha[4]={bondParam.phi0.x,
                         bondParam.phi0.y,
                         bondParam.phi0.z,
                         bondParam.phi0.w};

            real e=0;
            for(int n=1;n<=4;n++){

                Dihedral_::BondParameters bP;
                bP.n    = n;
                bP.K    = K[n-1];
                bP.phi0 = pha[n-1];

                e+=Dihedral_::energy(cos_dih,sin_dih,computational,bP);
            }

            return e;
        }

        static inline __device__ real energyDerivate(const real& cos_dih,const real& sin_dih,
                                                     const ComputationalData &computational,
                                                     const BondParameters &bondParam){

            real K[4]  ={bondParam.K.x,
                         bondParam.K.y,
                         bondParam.K.z,
                         bondParam.K.w};

            real pha[4]={bondParam.phi0.x,
                         bondParam.phi0.y,
                         bondParam.phi0.z,
                         bondParam.phi0.w};

            real fmod=0;
            for(int n=1;n<=4;n++){

                Dihedral_::BondParameters bP;
                bP.n    = n;
                bP.K    = K[n-1];
                bP.phi0 = pha[n-1];

                fmod+=Dihedral_::energyDerivate(cos_dih,sin_dih,computational,bP);
            }

            return fmod;
        }

      static inline __device__ real energySecondDerivate(const real& cos_dih,const real& sin_dih,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){

            real K[4]  ={bondParam.K.x,
                         bondParam.K.y,
                         bondParam.K.z,
                         bondParam.K.w};

            real pha[4]={bondParam.phi0.x,
                         bondParam.phi0.y,
                         bondParam.phi0.z,
                         bondParam.phi0.w};

            real derivfmod=0;
            for(int n=1;n<=4;n++){

                Dihedral_::BondParameters bP;
                bP.n    = n;
                bP.K    = K[n-1];
                bP.phi0 = pha[n-1];

                derivfmod+=Dihedral_::energySecondDerivate(cos_dih,sin_dih,computational,bP);
            }

            return derivfmod;
        }

    };

    using Dihedral                = Bond4Hessian_<Dihedral_>;
    using DihedralCommon_n_K_phi0 = Bond4Hessian_<DihedralCommon_n_K_phi0_>;
    using Dihedral4               = Bond4Hessian_<Dihedral4_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond4,Dihedral,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond4::Dihedral>
)

REGISTER_BOND_INTERACTOR(
    Bond4,DihedralCommon_n_K_phi0,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond4::DihedralCommon_n_K_phi0>
)

REGISTER_BOND_INTERACTOR(
    Bond4,Dihedral4,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond4::Dihedral4>
)
