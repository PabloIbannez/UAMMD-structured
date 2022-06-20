#ifndef __GAUSSIAN_BOND2__
#define __GAUSSIAN_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    struct Gaussian_{

        struct BondInfo{
            real E;
            real r0;
            real D;
        };

        Box box;

        struct Parameters{
            Box box;
        };

        Gaussian_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
        
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real E   = bi.E;
            const real r0  = bi.r0;
            const real D   = bi.D;

            const real r2 = dot(rij, rij);
            
            real3 f = CommonPotentials::GaussianWell::force(rij,r2,E,r0,D);

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
            
            const real E   = bi.E;
            const real r0  = bi.r0;
            const real D   = bi.D;
            
            const real r2 = dot(rij, rij);

            const real e = CommonPotentials::GaussianWell::energy(rij,r2,E,r0,D);
            
            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in>>bi.E>>bi.r0>>bi.D;
            return bi;
        }

    };
    
    struct GaussianConst_E_r0_D_ : public Gaussian_{
    
        struct BondInfo{};

        real E;
        real r0;
        real D;
        
        struct Parameters : Gaussian_::Parameters{
            real E;
            real r0;
            real D;
        };
        
        GaussianConst_E_r0_D_(Parameters param):Gaussian_(param),
                                                E(param.E),
                                                r0(param.r0),
                                                D(param.D){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Gaussian_::BondInfo biB = {E,r0,D};

            return Gaussian_::force(i,j,bond_index,posi,posj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Gaussian_::BondInfo biB = {E,r0,D};

            return Gaussian_::energy(i,j,bond_index,posi,posj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }
    };

    using Gaussian             = Bond2<addVirialStress<Gaussian_>>;
    using GaussianConst_E_r0_D = Bond2<addVirialStress<GaussianConst_E_r0_D_>>;

}}}}

#endif
