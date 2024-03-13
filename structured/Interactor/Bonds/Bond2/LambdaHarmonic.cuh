#ifndef __LAMBDA_HARMONIC_BOND2__
#define __LAMBDA_HARMONIC_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct LambdaHarmonic_{

        struct ComputationalData{
            real4* pos;
            Box    box;

            real lambda;
            int  n;
        };

        //Potential parameters

        struct StorageData{
            int  n;
        };

        struct BondParameters{

            real K;
            real r0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box = gd->getEnsemble()->getBox();

            computational.lambda = gd->getEnsemble()->getLambda();
            computational.n      = storage.n;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.n = data.getParameter<int>("n",2);

            System::log<System::MESSAGE>("[LambdaFixedHarmonic] n: %i",storage.n);

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.K    = bondParametersMap.at("K");
            param.r0   = bondParametersMap.at("r0");

            return param;
        }

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K   = bondParam.K;
            const real r0  = bondParam.r0;

            const real lambda = computational.lambda;
            const int  n      = computational.n;

            const real r2 = dot(rij, rij);

            const real e = BasicPotentials::lambdaHarmonic::energy(rij,r2,K,r0,lambda,n);

            return e;
        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K   = bondParam.K;
            const real r0  = bondParam.r0;

            const real lambda = computational.lambda;
            const int  n      = computational.n;

            const real r2 = dot(rij, rij);

            real3 f = BasicPotentials::lambdaHarmonic::force(rij,r2,K,r0,lambda,n);

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;
        }

        static inline __device__ real lambdaDerivative(int index_i, int index_j,
                                                       int currentParticleIndex,
                                                       const ComputationalData &computational,
                                                       const BondParameters   &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K   = bondParam.K;
            const real r0  = bondParam.r0;

            const real lambda = computational.lambda;
            const int  n      = computational.n;

            const real r2 = dot(rij, rij);

            const real ld = BasicPotentials::lambdaHarmonic::lambdaDerivative(rij,r2,K,r0,lambda,n);

            return ld;
        }
    };

    using LambdaHarmonic = Bond2Lambda_<LambdaHarmonic_>;

}}}}

#endif
