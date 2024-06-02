#ifndef __IDP_FOURIER_BOND4__
#define __IDP_FOURIER_BOND4__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond4{

    //Calibrated Langevin-dynamics simulations of intrinsically disordered proteins
    //W. Wendell Smith, Po-Yi Ho, and Corey S. O’Hern

    struct IDP_Fourier_{

        private:

            static constexpr real4 As_ = { 0.705, -0.313, -0.079, 0.041};
            static constexpr real4 Bs_ = {-0.175, -0.093,  0.030, 0.030};

        public:

            struct ComputationalData{
                real4* pos;
                Box    box;

                real4 As;
                real4 Bs;
            };

            struct StorageData{
                real4 As;
                real4 Bs;
            };

            struct BondParameters{};

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

                computational.As = storage.As;
                computational.Bs = storage.Bs;

                return computational;
            }

            //Storage data reader

            static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                       std::shared_ptr<ParticleGroup> pg,
                                                       DataEntry& data){

                StorageData storage;

                //Check if the used units are KcalMol_A
                std::string unitsName = gd->getUnits()->getSubType();
                if(unitsName != "KcalMol_A"){
                    System::log<System::CRITICAL>("[IDP_Fourier] The potential is only defined for KcalMol_A units."
                                                  " But the used units are %s", unitsName.c_str());
                }

                real kBT = gd->getUnits()->getBoltzmannConstant()*gd->getEnsemble()->getTemperature();

                storage.As = As_*kBT;
                storage.Bs = Bs_*kBT;

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

                const real As[4] = {computational.As.x, computational.As.y,
                                    computational.As.z, computational.As.w};

                const real Bs[4] = {computational.Bs.x, computational.Bs.y,
                                    computational.Bs.z, computational.Bs.w};

                real e = real(0.0);

                real cosnt;
                real sinnt;

                real tmp;

                for(int s = 1; s <= 4; ++s){

                    //cos((s+1)*phi) = cos(phi)*cos(s*phi) - sin(phi)*sin(s*phi)
                    //sin((s+1)*phi) = sin(phi)*cos(s*phi) + cos(phi)*sin(s*phi)

                    //s = 0
                    real cosst=real(1);
                    real sinst=real(0);

                    for(int i = 0; i < s; ++i){
                        tmp   = cosst*cos_dih - sinst*sin_dih;
                        sinst = sinst*cos_dih + cosst*sin_dih;
                        cosst = tmp;
                    }

                    // Now we have cosst = cos(s*theta) and sinst = sin(s*theta)
                    e += As[s-1]*cosst + Bs[s-1]*sinst;
                }


                return e;
            }

            static inline __device__ real energyDerivate(const real& cos_dih,const real& sin_dih,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){

                const real As[4] = {computational.As.x, computational.As.y,
                                    computational.As.z, computational.As.w};

                const real Bs[4] = {computational.Bs.x, computational.Bs.y,
                                    computational.Bs.z, computational.Bs.w};


                real de = real(0.0);

                real cosnt;
                real sinnt;

                real tmp;

                for(int s = 1; s <= 4; ++s){

                    //cos((s+1)*phi) = cos(phi)*cos(s*phi) - sin(phi)*sin(s*phi)
                    //sin((s+1)*phi) = sin(phi)*cos(s*phi) + cos(phi)*sin(s*phi)

                    //s = 0
                    real cosst=real(1);
                    real sinst=real(0);

                    for(int i = 0; i < s; ++i){
                        tmp   = cosst*cos_dih - sinst*sin_dih;
                        sinst = sinst*cos_dih + cosst*sin_dih;
                        cosst = tmp;
                    }

                    // Now we have cosst = cos(s*theta) and sinst = sin(s*theta)
                    // e += As[s-1]*cosst + Bs[s-1]*sinst;
                    de += -s*As[s-1]*sinst + s*Bs[s-1]*cosst;
                }

                return de;

            }

      static inline __device__ real energySecondDerivate(const real& cos_dih,const real& sin_dih,
							 const ComputationalData &computational,
							 const BondParameters &bondParam){

                const real As[4] = {computational.As.x, computational.As.y,
                                    computational.As.z, computational.As.w};

                const real Bs[4] = {computational.Bs.x, computational.Bs.y,
                                    computational.Bs.z, computational.Bs.w};


                real de = real(0.0);

                real cosnt;
                real sinnt;

                real tmp;

                for(int s = 1; s <= 4; ++s){

                    //cos((s+1)*phi) = cos(phi)*cos(s*phi) - sin(phi)*sin(s*phi)
                    //sin((s+1)*phi) = sin(phi)*cos(s*phi) + cos(phi)*sin(s*phi)

                    //s = 0
                    real cosst=real(1);
                    real sinst=real(0);

                    for(int i = 0; i < s; ++i){
                        tmp   = cosst*cos_dih - sinst*sin_dih;
                        sinst = sinst*cos_dih + cosst*sin_dih;
                        cosst = tmp;
                    }

                    // Now we have cosst = cos(s*theta) and sinst = sin(s*theta)
                    // e += As[s-1]*cosst + Bs[s-1]*sinst;
                    de += -s*s*As[s-1]*cosst - s*s*Bs[s-1]*sinst;
                }

                return de;

            }

    };

    using IDP_Fourier = Bond4Hessian_<IDP_Fourier_>;

}}}}

#endif
