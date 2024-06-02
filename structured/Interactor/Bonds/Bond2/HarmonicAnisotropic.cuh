#ifndef __HARMONIC_ANISOTROPIC_BOND2__
#define __HARMONIC_ANISOTROPIC_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct HarmonicAnisotropic_{

        struct ComputationalData{
            real4* pos;
            Box    box;
        };

        //Potential parameters

        struct StorageData{};

        struct BondParameters{

            real3 K;
            real3 r0;
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
            param.K    = bondParametersMap.at("K");
            param.r0   = bondParametersMap.at("r0");

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
	    real3 rij_abs = abs(rij);
	    
            const real3 K   = bondParam.K;
            const real3 r0  = bondParam.r0;

	    const real e = BasicPotentials::HarmonicAnisotropic::energy(rij,K,r0);
            return e;
        }


        static inline __device__ real3 force(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            real3 posi = make_real3(computational.pos[index_i]);
            real3 posj = make_real3(computational.pos[index_j]);

            const real3 rij = computational.box.apply_pbc(posj-posi);

            const real3 K   = bondParam.K;
            const real3 r0  = bondParam.r0;
	    real3 f = BasicPotentials::HarmonicAnisotropic::force(rij,K,r0);
	    
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

            const real3 K = bondParam.K;
	    tensor3 H = tensor3(0.0);


	    if (currentParticleIndex != index_i && currentParticleIndex != index_j){
		    H     = tensor3(0.0);

	    } else if        (currentParticleIndex == index_i){
	      H.xx = -K.x;
	      H.yy = -K.y;
	      H.zz = -K.z;
            } else if (currentParticleIndex == index_j){
	      H.xx = -K.x;
              H.yy = -K.y;
              H.zz = -K.z;
            }
	    
            return H;
        }

    };

    using HarmonicAnisotropic = Bond2Hessian_<HarmonicAnisotropic_>;
 
}}}}

#endif
