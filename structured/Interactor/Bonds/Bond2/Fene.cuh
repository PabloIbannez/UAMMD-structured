#ifndef __FENE_BOND2__
#define __FENE_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct Fene_{

        struct ComputationalData{
            real4* pos;
            Box    box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{

            real r0;
            real K;
            real R0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
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

            param.r0  = bondParametersMap.at("r0");
            param.K   = bondParametersMap.at("K");
            param.R0  = bondParametersMap.at("R0");

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

            const real K  = bondParam.K;
            const real r0 = bondParam.r0;
            const real R0 = bondParam.R0;

            const real r2 = dot(rij, rij);

            real e = BasicPotentials::Fene::energy(rij,r2,K,r0,R0);

            return e;
        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K  = bondParam.K;
            const real r0 = bondParam.r0;
            const real R0 = bondParam.R0;

            const real r2 = dot(rij, rij);

            real3 f = BasicPotentials::Fene::force(rij,r2,K,r0,R0);

            if        (currentParticleIndex == index_i){
            } else if (currentParticleIndex == index_j){
                f=-f;
            }

            return f;
        }

      static inline __device__ tensor3 hessian(int index_i, int index_j,
					                                     int currentParticleIndex,
					       const ComputationalData &computational,
					       const BondParameters &bondParam){

            tensor3 H = tensor3(0.0);

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real K  = bondParam.K;
            const real r0 = bondParam.r0;
            const real R0 = bondParam.R0;

            const real r2 = dot(rij, rij);

            if        (currentParticleIndex == index_i){
	      H = BasicPotentials::Fene::hessian(rij, r2, K, r0, R0);
            } else if (currentParticleIndex == index_j){
	      H = BasicPotentials::Fene::hessian(-rij, r2, K, r0, R0);
            }

            return H;
        }

    };

    struct FeneCommon_K_R0_{

        struct ComputationalData: public Fene_::ComputationalData{
            real K;
            real R0;
        };

        //Potential parameters

        struct StorageData: Fene_::StorageData{
            real K;
            real R0;
        };

        struct BondParameters{

            real r0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){
            ComputationalData computational;
            static_cast<Fene_::ComputationalData&>(computational) =
            Fene_::getComputationalData(gd,pg,storage,computables,st);

            computational.K  = storage.K;
            computational.R0 = storage.R0;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.K  = data.getParameter<real>("K");
            storage.R0 = data.getParameter<real>("R0");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.r0  = bondParametersMap.at("r0");

            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            Fene_::BondParameters bP;
            bP.R0 = computational.R0;
            bP.K  = computational.K;
            bP.r0 = bondParam.r0;

            return Fene_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            Fene_::BondParameters bP;
            bP.R0 = computational.R0;
            bP.K  = computational.K;
            bP.r0 = bondParam.r0;

            return Fene_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }


      static inline __device__ tensor3 hessian(int index_i, int index_j,
					       int currentParticleIndex,
					       const ComputationalData &computational,
					       const BondParameters &bondParam){


	Fene_::BondParameters bP;
	bP.R0 = computational.R0;
	bP.K  = computational.K;
	bP.r0 = bondParam.r0;

	return Fene_::hessian(index_i,index_j,currentParticleIndex,computational,bP);
        }
    };

    struct FeneCommon_r0_K_R0_{

        struct ComputationalData : public Fene_::ComputationalData {
            real K;
            real R0;
            real r0;
        };

        //Potential parameters

        struct StorageData: public Fene_::StorageData{
            real K;
            real R0;
            real r0;
        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;
            static_cast<Fene_::ComputationalData&>(computational) =
            Fene_::getComputationalData(gd,pg,storage,computables,st);

            computational.K  = storage.K;
            computational.R0 = storage.R0;
            computational.r0 = storage.r0;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.K  = data.getParameter<real>("K");
            storage.R0 = data.getParameter<real>("R0");
            storage.r0 = data.getParameter<real>("r0");

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
                                             const BondParameters   &bondParam){

            Fene_::BondParameters bP;
            bP.R0 = computational.R0;
            bP.K  = computational.K;
            bP.r0 = computational.r0;

            return Fene_::energy(index_i,index_j,currentParticleIndex,computational,bP);
        }

        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters   &bondParam){

            Fene_::BondParameters bP;
            bP.R0 = computational.R0;
            bP.K  = computational.K;
            bP.r0 = computational.r0;

            return Fene_::force(index_i,index_j,currentParticleIndex,computational,bP);

        }

      static inline __device__ tensor3 hessian(int index_i, int index_j,
              int currentParticleIndex,
              const ComputationalData &computational,
              const BondParameters &bondParam){


          Fene_::BondParameters bP;
          bP.R0 = computational.R0;
          bP.K  = computational.K;
          bP.r0 = computational.r0;

          return Fene_::hessian(index_i,index_j,currentParticleIndex,computational,bP);
      }
    };

    using Fene               = Bond2Hessian_<Fene_>;
    using FeneCommon_K_R0    = Bond2Hessian_<FeneCommon_K_R0_>;
    using FeneCommon_r0_K_R0 = Bond2Hessian_<FeneCommon_r0_K_R0_>;

}}}}

#endif
