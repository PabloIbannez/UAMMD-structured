#ifndef __FENE_BOND2__
#define __FENE_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    struct Fene_{

        struct BondInfo{
            real K;
            real R02;
            real r0;
        };
        
        Box box;
        
        struct Parameters{
            Box box;
        };

        Fene_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
        
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real K   = bi.K;
            const real R02 = bi.R02;
            const real r0  = bi.r0;

            const real  r  = sqrt(dot(rij, rij));
            const real dr  = r-r0;
            const real dr2 = dr*dr;

            const real fmod = K*dr/(real(1.0)-dr2/R02);

            real3 f = fmod*rij/r;
            
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
            
            real3 rij = box.apply_pbc(posj-posi);
            
            const real K   = bi.K;
            const real R02 = bi.R02;
            const real r0  = bi.r0;
            
            real  r  = sqrt(dot(rij, rij));
            real dr  = r-r0;
            real dr2 = dr*dr;

            real e = -(K*R02/real(2.0))*log(real(1.0)-dr2/R02);
            
            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in>>bi.K>>bi.R02>>bi.r0;
            bi.R02=bi.R02*bi.R02;
            return bi;
        }
    };

    struct FeneConst_K_R0_ : public Fene_{
        
        struct BondInfo{
            real r0;
        };
        
        real K;
        real R02;

        struct Parameters : Fene_::Parameters{
            real K;
            real R0;
        };

        FeneConst_K_R0_(Parameters par):Fene_(par),
                                        K(par.K),
                                        R02(par.R0*par.R0){}
        
        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Fene_::BondInfo biB = {K,R02,bi.r0};

            return Fene_::force(i,j,bond_index,posi,posj,biB);
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi) {            
            
            const Fene_::BondInfo biB = {K,R02,bi.r0};

            return Fene_::energy(i,j,bond_index,posi,posj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in>>bi.r0;
            return bi;
        }

    };

    using Fene           = Bond2<addVirial<Fene_>>;
    using FeneConst_K_R0 = Bond2<addVirial<FeneConst_K_R0_>>;

}}}}

#endif
