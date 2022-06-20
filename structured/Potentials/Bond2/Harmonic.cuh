#ifndef __HARMONIC_BOND2__
#define __HARMONIC_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    struct Harmonic_{

        struct BondInfo{
            real K;
            real r0;
        };

        Box box;

        struct Parameters{
            Box box;
        };

        Harmonic_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
        
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real K   = bi.K;
            const real r0  = bi.r0;

            const real r2 = dot(rij, rij);
            
            real3 f = CommonPotentials::Harmonic::force(rij,r2,K,r0);

            if        (bond_index == i){
                /*
                if(sqrt(r2)>r0+5.0 or sqrt(r2)<r0-5.0){
                    printf("i:%i j:%i ,r:%f r0:%f, pi: %f %f %f, pj: %f %f %f\n",i,j,
                                                                                 sqrt(r2),r0,
                                                                                 posi.x,posi.y,posi.z, 
                                                                                 posj.x,posj.y,posj.z);
                }*/
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
            
            const real K   = bi.K;
            const real r0  = bi.r0;
            
            const real r2 = dot(rij, rij);

            const real e = CommonPotentials::Harmonic::energy(rij,r2,K,r0);
            
            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in>>bi.r0>>bi.K;
            return bi;
        }

    };
    
    struct HarmonicConst_K_ : public Harmonic_{
    
        struct BondInfo{
            real r0;
        };

        real K;
        
        struct Parameters : Harmonic_::Parameters{
            real K;
        };
        
        HarmonicConst_K_(Parameters param):Harmonic_(param),
                                           K(param.K){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Harmonic_::BondInfo biB = {K,bi.r0};

            return Harmonic_::force(i,j,bond_index,posi,posj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Harmonic_::BondInfo biB = {K,bi.r0};

            return Harmonic_::energy(i,j,bond_index,posi,posj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in>>bi.r0;
            return bi;
        }
    };

    struct HarmonicConst_K_r0_ : public Harmonic_{
    
        struct BondInfo{};

        real K;
        real r0;
        
        struct Parameters : Harmonic_::Parameters{
            real K;
            real r0;
        };
        
        HarmonicConst_K_r0_(Parameters param):Harmonic_(param),
                                              K(param.K),
                                              r0(param.r0){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Harmonic_::BondInfo biB = {K,r0};

            return Harmonic_::force(i,j,bond_index,posi,posj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Harmonic_::BondInfo biB = {K,r0};

            return Harmonic_::energy(i,j,bond_index,posi,posj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }
    };

    using Harmonic           = Bond2<addVirialStress<Harmonic_>>;
    using HarmonicConst_K    = Bond2<addVirialStress<HarmonicConst_K_>>;
    using HarmonicConst_K_r0 = Bond2<addVirialStress<HarmonicConst_K_r0_>>;

}}}}

#endif
