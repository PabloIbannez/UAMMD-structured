#ifndef __FIXED_HARMONIC_BOND1__
#define __FIXED_HARMONIC_BOND1__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond1{

    struct FixedHarmonic_{

        struct ComputationalData{
            real4* pos;
            Box    box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{
            real K;
            real r0;
            real3 pos;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage){

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
	    const real r2    = dot(rij, rij);
	    
            const real K   = bondParam.K;
            const real r0  = bondParam.r0;

            return BasicPotentials::Harmonic::force(rij,r2,K,r0);
        }

        static inline __device__ real energy(int index_i,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

	  const real3 posi = make_real3(computational.pos[index_i]);
            const real3 posj = bondParam.pos;
            const real3 rij  = computational.box.apply_pbc(posj-posi);
	    const real r2    = dot(rij, rij);
	    
            const real K  = bondParam.K;
            const real r0 = bondParam.r0;

            const real e = BasicPotentials::Harmonic::energy(rij,r2,K,r0);

            return e;
        }

      static inline __device__ tensor3 hessian(int index_i,
					       int currentParticleIndex,
					       const ComputationalData &computational,
					       const BondParameters &bondParam){
	
	tensor3 H = tensor3();
	const real3 posi = make_real3(computational.pos[index_i]);
	const real3 posj = bondParam.pos;
	const real3 rij  = computational.box.apply_pbc(posj-posi);
	const real r2    = dot(rij, rij);
	
	const real K   = bondParam.K;
	const real r0  = bondParam.r0;
	
        H = -BasicPotentials::Harmonic::hessian(rij,r2,K,r0);
	return H;
      }

    };

    using FixedHarmonic = Bond1Hessian_<FixedHarmonic_>;

}}}}

#endif
