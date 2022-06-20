#ifndef __HARMONIC_ANGULAR__
#define __HARMONIC_ANGULAR__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond3{
            
    struct HarmonicAngular_{
    
        private:
            
            Box box;

        public:

            struct Parameters{
                Box box;
            };

            HarmonicAngular_(Parameters par): box(par.box){}

            struct BondInfo{
                real ang0;
                real K;
            };

            struct staticFunctions{
                static inline __device__ real f(const real& ang,const BondInfo& bi){
                    real adiff = ang - bi.ang0;
                    return -bi.K*adiff;
                }
                
                static inline __device__ real e(const real& ang,const BondInfo& bi){
                    real adiff = ang - bi.ang0;
                    return real(0.5)*bi.K*adiff*adiff; 
                }
            };

            inline __device__ real3 force(int i, int j, int k,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const BondInfo &bi){

                real3 fi;
                real3 fk;

                angularForce<BondInfo,staticFunctions>(posi,posj,posk,
                                                       box,
                                                       bi,
                                                       fi,fk);

                if(bond_index==i){
                    return -fi;
                }
                else if(bond_index==j){
                    return fi+fk;
                }
                else if(bond_index==k){
                    return -fk;
                }
                
                return make_real3(0);
            }
            
            inline __device__ real virial(int i, int j, int k,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const BondInfo &bi){

                return real(0);
            }
            
            inline __device__ tensor3 stress(int i, int j, int k,
                                             int bond_index,
                                             const real3 &posi,
                                             const real3 &posj,
                                             const real3 &posk,
                                             const BondInfo &bi){

                return tensor3(0);
            }

            inline __device__ real energy(int i, int j, int k,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const BondInfo &bi){

                real e;

                angularEnergy<BondInfo,staticFunctions>(posi,posj,posk,
                                                        box,
                                                        bi,
                                                        e);
                
                return e/real(3.0);
            }

            static BondInfo readBond(std::istream &in){
                BondInfo bi;
                in>>bi.ang0>>bi.K;
                return bi;
            }

    };
    
    struct HarmonicAngularConst_K_ : public HarmonicAngular_{
    
            struct BondInfo{
                real ang0;
            };

            real K;

            struct Parameters : public HarmonicAngular_::Parameters {
                real K;
            };

            HarmonicAngularConst_K_(Parameters par):HarmonicAngular_(par),
                                                    K(par.K){}

            inline __device__ real3 force(int i, int j, int k,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const BondInfo &bi){

                const HarmonicAngular_::BondInfo biB = {.ang0=bi.ang0,
                                                        .K=K,};
            
                return HarmonicAngular_::force(i,j,k,
                                               bond_index,
                                               posi,
                                               posj,
                                               posk,
                                               biB);
                
                return make_real3(0);
            }
            
            inline __device__ tensor3 virial(int i, int j, int k,
                                             int bond_index,
                                             const real3 &posi,
                                             const real3 &posj,
                                             const real3 &posk,
                                             const BondInfo &bi){

                return tensor3(0);
            }
            
            inline __device__ tensor3 stress(int i, int j, int k,
                                             int bond_index,
                                             const real3 &posi,
                                             const real3 &posj,
                                             const real3 &posk,
                                             const BondInfo &bi){

                return tensor3(0);
            }

            inline __device__ real energy(int i, int j, int k,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const BondInfo &bi){

                const HarmonicAngular_::BondInfo biB = {.ang0=bi.ang0,
                                                        .K=K,};
            
                return HarmonicAngular_::energy(i,j,k,
                                                bond_index,
                                                posi,
                                                posj,
                                                posk,
                                                biB);
            }

            static BondInfo readBond(std::istream &in){
                BondInfo bi;
                in >> bi.ang0;
                return bi;
            }

    };
    
    struct HarmonicAngularConst_K_ang0_ : public HarmonicAngular_{
    
            struct BondInfo{};

            real ang0;
            real K;

            struct Parameters : public HarmonicAngular_::Parameters {
                real ang0;
                real K;
            };

            HarmonicAngularConst_K_ang0_(Parameters par):HarmonicAngular_(par),
                                                         ang0(par.ang0),
                                                         K(par.K){}

            inline __device__ real3 force(int i, int j, int k,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const BondInfo &bi){

                const HarmonicAngular_::BondInfo biB = {.ang0=ang0,
                                                        .K=K,};
            
                return HarmonicAngular_::force(i,j,k,
                                               bond_index,
                                               posi,
                                               posj,
                                               posk,
                                               biB);
                
                return make_real3(0);
            }
            
            inline __device__ tensor3 virial(int i, int j, int k,
                                             int bond_index,
                                             const real3 &posi,
                                             const real3 &posj,
                                             const real3 &posk,
                                             const BondInfo &bi){

                return tensor3(0);
            }
            
            inline __device__ tensor3 stress(int i, int j, int k,
                                             int bond_index,
                                             const real3 &posi,
                                             const real3 &posj,
                                             const real3 &posk,
                                             const BondInfo &bi){

                return tensor3(0);
            }

            inline __device__ real energy(int i, int j, int k,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const BondInfo &bi){

                const HarmonicAngular_::BondInfo biB = {.ang0=ang0,
                                                        .K=K,};
            
                return HarmonicAngular_::energy(i,j,k,
                                                bond_index,
                                                posi,
                                                posj,
                                                posk,
                                                biB);
            }

            static BondInfo readBond(std::istream &in){
                BondInfo bi;
                return bi;
            }

    };

    using HarmonicAngular             = Bond3<HarmonicAngular_>;
    using HarmonicAngularConst_K      = Bond3<HarmonicAngularConst_K_>;
    using HarmonicAngularConst_K_ang0 = Bond3<HarmonicAngularConst_K_ang0_>;

}}}}

#endif
