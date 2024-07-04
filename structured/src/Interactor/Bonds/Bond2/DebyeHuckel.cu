#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond2/Bond2.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct DebyeHuckel_{

        struct ComputationalData{
            real4* pos;
            Box box;
            real ELECOEF;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{

            real chgProduct;
            real dielectricConstant;
            real debyeLength;
            real cutOff;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos     = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box     = gd->getEnsemble()->getBox();
            computational.ELECOEF = gd->getUnits()->getElectricConversionFactor();

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

            param.chgProduct         = bondParametersMap.at("chgProduct");
            param.dielectricConstant = bondParametersMap.at("dielectricConstant");
            param.debyeLength        = bondParametersMap.at("debyeLength");
            param.cutOff             = bondParametersMap.at("cutOff");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            real3 rij = computational.box.apply_pbc(posj-posi);

            real  r2 = dot(rij, rij);

            real e=real(0.0);
            if(r2<(bondParam.cutOff*bondParam.cutOff)){
                e = BasicPotentials::DebyeHuckel::DebyeHuckel::energy(rij,r2,
                                                                      computational.ELECOEF,
                                                                      bondParam.chgProduct,
                                                                      bondParam.dielectricConstant,
                                                                      bondParam.debyeLength);
            }

            return e;
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            real3 rij = computational.box.apply_pbc(posj-posi);

            real  r2 = dot(rij, rij);

            real3 f=make_real3(0.0);
            if(r2<(bondParam.cutOff*bondParam.cutOff)){
                f = BasicPotentials::DebyeHuckel::DebyeHuckel::force(rij,r2,
                                                                     computational.ELECOEF,
                                                                     bondParam.chgProduct,
                                                                     bondParam.dielectricConstant,
                                                                     bondParam.debyeLength);

                if        (currentParticleIndex == index_i){
                } else if (currentParticleIndex == index_j){
                    f=-f;
                }
            }

            return f;

        }


    };

    using DebyeHuckel = Bond2_<DebyeHuckel_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond2,DebyeHuckel,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::DebyeHuckel>
)
