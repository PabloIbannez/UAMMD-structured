#ifndef __HARMONIC_ANGULAR__
#define __HARMONIC_ANGULAR__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond3{

    struct HarmonicAngular_{

        //Potential parameters

        struct ComputationalData{
            Box    box;
            real4* pos;
        };

        struct StorageData{};

        struct BondParameters{
            real K;
            real theta0;
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.box = gd->getEnsemble()->getBox();
            computational.pos = pd->getPos(access::location::gpu,
                                           access::mode::read).raw();

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

            param.K      = bondParametersMap.at("K");
            param.theta0 = bondParametersMap.at("theta0");

            return param;
        }

      //Energy and force definition
      
      static inline __device__ real energy(const real& ang,
					   const ComputationalData &computational,
					   const BondParameters &bondParam){
	
	real adiff = ang - bondParam.theta0;
	return real(0.5)*bondParam.K*adiff*adiff;
      }

      static inline __device__ real energyDerivate(const real& ang,
						   const ComputationalData &computational,
						   const BondParameters &bondParam){
	
	real adiff = ang - bondParam.theta0;
	return bondParam.K*adiff;
      }

      static inline __device__ real energySecondDerivate(const real& ang,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){
	
	return bondParam.K;
      }

      
    };

    struct HarmonicAngularCommon_K_{

        //Potential parameters

        struct ComputationalData: public HarmonicAngular_::ComputationalData{
            real K;
        };

        struct StorageData: public HarmonicAngular_::StorageData{
            real K;
        };

        struct BondParameters{
            real theta0;
        };

        //Storage data reader

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage){
            ComputationalData computational;
            static_cast<HarmonicAngular_::ComputationalData&>(computational) = HarmonicAngular_::getComputationalData(gd, pg, storage);

            computational.K = storage.K;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.K   = data.getParameter<real>("K");

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;

            param.theta0 = bondParametersMap.at("theta0");

            return param;
        }

        //Angular potential functions

        static inline __device__ real energy(const real& ang,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            HarmonicAngular_::BondParameters bP;
            bP.theta0 = bondParam.theta0;
            bP.K      = computational.K;

            return HarmonicAngular_::energy(ang,computational,bP);
        }

        static inline __device__ real energyDerivate(const real& ang,
                                                     const ComputationalData &computational,
                                                     const BondParameters &bondParam){

            HarmonicAngular_::BondParameters bP;
            bP.theta0 = bondParam.theta0;
            bP.K      = computational.K;

            return HarmonicAngular_::energyDerivate(ang,computational,bP);
        }

      static inline __device__ real energySecondDerivate(const real& ang,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){
	return computational.K;
      }
    };
  
    struct HarmonicAngularCommon_K_theta0_{

        //Potential parameters

        struct ComputationalData: public HarmonicAngular_::ComputationalData{
            real K;
            real theta0;
        };

        struct StorageData: public HarmonicAngular_::StorageData{
            real K;
            real theta0;
        };

        struct BondParameters{};

        //Storage data reader

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage){
            ComputationalData computational;
            static_cast<HarmonicAngular_::ComputationalData&>(computational) =
            HarmonicAngular_::getComputationalData(gd, pg, storage);

            computational.K      = storage.K;
            computational.theta0 = storage.theta0;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.K      = data.getParameter<real>("K");
            storage.theta0 = data.getParameter<real>("theta0");

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

        static inline __device__ real energy(const real& ang,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            HarmonicAngular_::BondParameters bP;
            bP.theta0 = computational.theta0;
            bP.K      = computational.K;

            return HarmonicAngular_::energy(ang,computational,bP);
        }

      static inline __device__ real energyDerivate(const real& ang,
                                                     const ComputationalData &computational,
                                                     const BondParameters &bondParam){

            HarmonicAngular_::BondParameters bP;
            bP.theta0 = computational.theta0;
            bP.K      = computational.K;

            return HarmonicAngular_::energyDerivate(ang,computational,bP);
        }

      static inline __device__ real energySecondDerivate(const real& ang,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){
	
	return computational.K;
      }
    };

    using HarmonicAngular                = Bond3Hessian_<HarmonicAngular_>;
    using HarmonicAngularCommon_K        = Bond3Hessian_<HarmonicAngularCommon_K_>;
    using HarmonicAngularCommon_K_theta0 = Bond3Hessian_<HarmonicAngularCommon_K_theta0_>;

}}}}

#endif
