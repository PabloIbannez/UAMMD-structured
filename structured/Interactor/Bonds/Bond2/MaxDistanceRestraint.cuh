#ifndef __MAX_DISTANCE_RESTRAINT_BOND2__
#define __MAX_DISTANCE_RESTRAINT_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct MaxDistanceRestraint_{

        struct ComputationalData{
            real4* pos;
            Box box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{
            real maxDistance;
            real K;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage){
            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box = gd->getEnsemble()->getBox(); //We ask global data for the box

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

            //Example:
            param.K           = bondParametersMap.at("K");
            param.maxDistance = bondParametersMap.at("maxDistance");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K           = bondParam.K;
            const real maxDistance = bondParam.maxDistance;

            const real r2 = dot(rij, rij);

            real e = 0;
            if(r2 > maxDistance*maxDistance){
                e = BasicPotentials::Harmonic::energy(rij,r2,K,maxDistance);
            }

            return e;
        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K           = bondParam.K;
            const real maxDistance = bondParam.maxDistance;

            const real r2 = dot(rij, rij);

            real3 f = make_real3(0);

            if(r2 > maxDistance*maxDistance){
                f = BasicPotentials::Harmonic::force(rij,r2,K,maxDistance);
            }

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;
        }

    };

    using MaxDistanceRestraint = Bond2_<MaxDistanceRestraint_>;

}}}}

#endif
