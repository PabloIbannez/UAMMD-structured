#ifndef __EXPONENTIAL_ANGULAR__
#define __EXPONENTIAL_ANGULAR__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond3{
            
    struct BestChenHummerAngular_{
    
        private:
            
            static constexpr real gamma = real(0.1);

            static constexpr real k_alpha = real(106.4);
            static constexpr real k_beta  = real(26.3);
            
            static constexpr real theta_alpha = real(1.60);
            static constexpr real theta_beta  = real(2.27);

            static constexpr real epsilon_alpha = real(4.3);

            Box box;

        public:

            struct Parameters{
                Box box;
            };

            BestChenHummerAngular_(Parameters par):box(par.box){}

            struct BondInfo{};

            struct staticFunctions{
                static inline __device__ real f(const real& ang,const BondInfo& bi){
                    
                    const real adiff_alpha  = ang-theta_alpha;
                    const real adiff_alpha2 = adiff_alpha*adiff_alpha;
                    const real exp_alpha = exp(-gamma*(k_alpha*adiff_alpha2+epsilon_alpha));
                    
                    const real adiff_beta  = ang-theta_beta;
                    const real adiff_beta2 = adiff_beta*adiff_beta;
                    const real exp_beta = exp(-gamma*k_beta*adiff_beta2);
                    
                    return -real(2.0)*(k_alpha*adiff_alpha*exp_alpha+k_beta*adiff_beta*exp_beta)/(exp_alpha+exp_beta);
                }
                
                static inline __device__ real e(const real& ang,const BondInfo& bi){
                    
                    real adiff_alpha2 = ang-theta_alpha;
                         adiff_alpha2 = adiff_alpha2*adiff_alpha2;
                    const real exp_alpha = exp(-gamma*(k_alpha*adiff_alpha2+epsilon_alpha));
                    
                    real adiff_beta2 = ang-theta_beta;
                         adiff_beta2 = adiff_beta2*adiff_beta2;
                    const real exp_beta = exp(-gamma*k_beta*adiff_beta2);
                    
                    return -real(1.0/gamma)*log(exp_alpha+exp_beta); 
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

                real3 fi;
                real3 fk;

                angularForce<BondInfo,staticFunctions>(posi,posj,posk,
                                                       box,
                                                       bi,
                                                       fi,fk);
                
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
                return bi;
            }

    };

    using BestChenHummerAngular = Bond3<BestChenHummerAngular_>;

}}}}

#endif
