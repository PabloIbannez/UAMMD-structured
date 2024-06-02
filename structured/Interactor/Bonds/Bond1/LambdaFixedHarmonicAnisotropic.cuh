#ifndef __LAMBDA_FIXED_HARMONIC_ANISOTROPIC_BOND1__
#define __LAMBDA_FIXED_HARMONIC_ANISOTROPIC_BOND1__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond1{

    struct LambdaFixedHarmonicAnisotropic_{

        struct ComputationalData{
            real4* pos;
            Box    box;

            real lambda;
            int  n;
        };

        //Potential parameters

        struct StorageData{
            int n;
        };

        struct BondParameters{
            real3 K;
            real3 r0;
            real3 pos;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.box    = gd->getEnsemble()->getBox();
            computational.lambda = gd->getEnsemble()->getLambda();

            computational.n      = storage.n;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.n = data.getParameter<int>("n",2);

            System::log<System::MESSAGE>("[LambdaFixedHarmonicAnisotropic] n: %i",storage.n);

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.K    = bondParametersMap.at("K");
            param.r0   = bondParametersMap.at("r0");
            param.pos  = bondParametersMap.at("position");

            return param;
        }

        //Energy and force definition

        static inline __device__ real3 force(int index_i,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            const real3 posi = make_real3(computational.pos[index_i]);
            const real3 posj = bondParam.pos;
            const real3 rij  = computational.box.apply_pbc(posj-posi);

            const real3 K   = bondParam.K;
            const real3 r0  = bondParam.r0;

            const real lambda = computational.lambda;
            const int  n      = computational.n;

            return BasicPotentials::lambdaHarmonicAnisotropic::force(rij,K,r0,lambda,n);
        }

        static inline __device__ real energy(int index_i,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            const real3 posi = make_real3(computational.pos[index_i]);
            const real3 posj = bondParam.pos;
            const real3 rij  = computational.box.apply_pbc(posj-posi);

            const real3 K   = bondParam.K;
            const real3 r0  = bondParam.r0;

            const real lambda = computational.lambda;
            const int  n      = computational.n;

            const real e = BasicPotentials::lambdaHarmonicAnisotropic::energy(rij,K,r0,lambda,n);

            return e;
        }

        static inline __device__ real lambdaDerivative(int index_i,
                                                       int currentParticleIndex,
                                                       const ComputationalData &computational,
                                                       const BondParameters &bondParam){

            const real3 posi = make_real3(computational.pos[index_i]);
            const real3 posj = bondParam.pos;
            const real3 rij  = computational.box.apply_pbc(posj-posi);

            const real3 K   = bondParam.K;
            const real3 r0  = bondParam.r0;

            const real lambda = computational.lambda;
            const int  n      = computational.n;

            const real ld = BasicPotentials::lambdaHarmonicAnisotropic::lambdaDerivative(rij,K,r0,lambda,n);

            return ld;
        }

    };

    using LambdaFixedHarmonicAnisotropic = Bond1Lambda_<LambdaFixedHarmonicAnisotropic_>;

}}}}

#endif
