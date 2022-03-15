#ifndef __KRATKY_POROD_ANGLE__
#define __KRATKY_POROD_ANGLE__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond3{
            
    struct KratkyPorod_{
    
        private:
            
            Box box;

        public:

            struct Parameters{
                Box box;
            };

            KratkyPorod_(Parameters par): box(par.box){}

            struct BondInfo{
                real K;
            };

            struct staticFunctions{
                static inline __device__ real f(const real& ang,const BondInfo& bi){
                    return bi.K*sin(ang);
                }
                
                static inline __device__ real e(const real& ang,const BondInfo& bi){
                    return bi.K*(real(1.0)+cos(ang));
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
            
            inline __device__ tensor3 virial(int i, int j, int k,
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
                in>>bi.K;
                return bi;
            }

    };
    
    struct KratkyPorodConst_K_ : public KratkyPorod_{
    
        struct BondInfo{};
        
        real K;

        struct Parameters: KratkyPorod_::Parameters{
            real K;
        };

        KratkyPorodConst_K_(Parameters par): KratkyPorod_(par),
                                             K(par.K){}

        inline __device__ real3 force(int i, int j, int k,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const real3 &posk,
                                      const BondInfo &bi){


            const KratkyPorod_::BondInfo biB = {.K=K};

            return KratkyPorod_::force(i,j,k,bond_index,posi,posj,posk,biB);

        }
        
        inline __device__ tensor3 virial(int i, int j, int k,
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

            const KratkyPorod_::BondInfo biB = {.K=K};

            return KratkyPorod_::energy(i,j,k,bond_index,posi,posj,posk,biB);
        }

        static BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }

    };

    using KratkyPorod        = Bond3<KratkyPorod_>;
    using KratkyPorodConst_K = Bond3<KratkyPorodConst_K_>;

}}}}

#endif
