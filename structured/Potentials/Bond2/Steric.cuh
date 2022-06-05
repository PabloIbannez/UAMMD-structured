#ifndef __STERIC_BOND2__
#define __STERIC_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    template<int power>
    struct Steric_{

        struct BondInfo{
            real epsilon;
            real sigma;
        };

        Box box;
        real cutOff2;

        struct Parameters{
            Box box;
            real cutOff;
        };

        Steric_(Parameters par):cutOff2(par.cutOff*par.cutOff),
                                box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
        
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            
            const real r2 = dot(rij, rij);

            real3 f = make_real3(0);

            if(r2<=cutOff2){
                f = CommonPotentials::Steric::Steric::force<power>(rij,r2,epsilon,sigma);
            } 
            
            if        (bond_index == i){
            } else if (bond_index == j){
                f=-f;
            }

            return f;
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi) {            
            
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            
            const real r2 = dot(rij, rij);

            real e = real(0);

            if(r2<=cutOff2){
                e = CommonPotentials::Steric::Steric::energy<power>(rij,r2,epsilon,sigma);
            } 
            
            return e;
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in>>bi.epsilon>>bi.sigma;
            return bi;
        }

    };

    template<int power>
    struct StericConst_epsilon_sigma_ : public Steric_<power>{
    
        struct BondInfo{};

        real epsilon;
        real sigma;
        
        struct Parameters : Steric_<power>::Parameters{
            real epsilon;
            real sigma;
        };
        
        StericConst_epsilon_sigma_(Parameters param):Steric_<power>(param),
                                                     epsilon(param.epsilon),
                                                     sigma(param.sigma){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const typename Steric_<power>::BondInfo biB = {epsilon,sigma};

            return Steric_<power>::force(i,j,bond_index,posi,posj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const typename Steric_<power>::BondInfo biB = {epsilon,sigma};

            return Steric_<power>::energy(i,j,bond_index,posi,posj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }
    };

    using Steric6                    = Bond2<addVirial<Steric_<6>>>;
    using Steric6Const_epsilon_sigma = Bond2<addVirial<StericConst_epsilon_sigma_<6>>>;
    
    using Steric12                    = Bond2<addVirial<Steric_<12>>>;
    using Steric12Const_epsilon_sigma = Bond2<addVirial<StericConst_epsilon_sigma_<12>>>;

}}}}

#endif
