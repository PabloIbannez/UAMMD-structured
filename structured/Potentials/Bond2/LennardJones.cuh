#ifndef __LENNARD_JONES_BOND2__
#define __LENNARD_JONES_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    template<class LennardJonesType>
    struct LennardJones_{

        struct BondInfo{
            real epsilon;
            real sigma;
        };
        
        Box box;

        struct Parameters{
            Box box;
        };

        LennardJones_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
    
            real3 rij = box.apply_pbc(posj-posi);
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            
            real  r2 = dot(rij, rij);
            
            real3 f = LennardJonesType::force(rij,r2,epsilon,sigma);
            
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
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            
            real  r2 = dot(rij, rij);
            
            real e = LennardJonesType::energy(rij,r2,epsilon,sigma);

            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            
            BondInfo bi;
            
            real sigma;
            real epsilon;
            in>>sigma>>epsilon;
            
            bi.sigma   = sigma;
            bi.epsilon = epsilon;
            
            return bi;
        }

    };
    
    template<class LennardJonesType>
    struct LennardJonesConst_e_ : public LennardJones_<LennardJonesType>{

        struct BondInfo{
            real sigma;
        };

        real epsilon;

        struct Parameters : LennardJones_<LennardJonesType>::Parameters{
            real epsilon;
        };

        LennardJonesConst_e_(Parameters par):LennardJones_<LennardJonesType>(par),
                                             epsilon(par.epsilon){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
    
            const typename LennardJones_<LennardJonesType>::BondInfo biB = {epsilon,bi.sigma};

            return LennardJones_<LennardJonesType>::force(i,j,bond_index,posi,posj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi) {
            
            const typename LennardJones_<LennardJonesType>::BondInfo biB = {epsilon,bi.sigma};

            return LennardJones_<LennardJonesType>::energy(i,j,bond_index,posi,posj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            
            BondInfo bi;
            
            real sigma;
            in>>sigma;
            
            bi.sigma   = sigma;
            
            return bi;
        }

    };
    
    struct LennardJonesKaranicolasBrooks_{

        struct BondInfo{
            real epsilon;
            real sigma;
        };
        
        Box box;

        struct Parameters{
            Box box;
        };

        LennardJonesKaranicolasBrooks_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
    
            real3 rij = box.apply_pbc(posj-posi);
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            
            real  r2 = dot(rij, rij);
            
            real3 f = CommonPotentials::LennardJones::KaranicolasBrooks::force(rij,r2,epsilon,sigma);
            
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
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            
            real  r2 = dot(rij, rij);
            
            real e = CommonPotentials::LennardJones::KaranicolasBrooks::energy(rij,r2,epsilon,sigma);

            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            
            BondInfo bi;
            
            real sigma;
            real epsilon;
            in>>sigma>>epsilon;
            
            bi.sigma   = sigma;
            bi.epsilon = epsilon;
            
            return bi;
        }

    };
    
    struct LennardJonesGaussian_{

        struct BondInfo{
            real epsilon;
            real sigma;
            real D;
        };
        
        Box box;

        struct Parameters{
            Box box;
        };

        LennardJonesGaussian_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
    
            real3 rij = box.apply_pbc(posj-posi);
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            const real D       = bi.D;
            
            real  r2 = dot(rij, rij);
            
            real3 f = CommonPotentials::ModifiedLennardJones::Gaussian::force(rij,r2,epsilon,sigma,D);
            
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
            
            const real epsilon = bi.epsilon;
            const real sigma   = bi.sigma;
            const real D       = bi.D;
            
            real r2 = dot(rij, rij);
            
            real e = CommonPotentials::ModifiedLennardJones::Gaussian::energy(rij,r2,epsilon,sigma,D);

            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            
            BondInfo bi;
            
            real sigma;
            real epsilon;
            real D;
            in>>sigma>>epsilon>>D;
            
            bi.sigma   = sigma;
            bi.epsilon = epsilon;
            bi.D       = D;
            
            return bi;
        }

    };
    
    struct LennardJonesGaussianConst_e_D_ : public LennardJonesGaussian_{

        struct BondInfo{
            real sigma;
        };

        real epsilon;
        real D;

        struct Parameters : LennardJonesGaussian_::Parameters{
            real epsilon;
            real D;
        };

        LennardJonesGaussianConst_e_D_(Parameters par):LennardJonesGaussian_(par),
                                                       epsilon(par.epsilon),
                                                       D(par.D){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
    
            const typename LennardJonesGaussian_::BondInfo biB = {epsilon,bi.sigma,D};

            return LennardJonesGaussian_::force(i,j,bond_index,posi,posj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi) {
            
            const typename LennardJonesGaussian_::BondInfo biB = {epsilon,bi.sigma,D};

            return LennardJonesGaussian_::energy(i,j,bond_index,posi,posj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            
            BondInfo bi;
            
            real sigma;
            in>>sigma;
            
            bi.sigma   = sigma;
            
            return bi;
        }

    };
        
    using LennardJonesType2        = Bond2<addVirial<LennardJones_<CommonPotentials::LennardJones::Type2>>>;
    using LennardJonesType2Const_e = Bond2<addVirial<LennardJonesConst_e_<CommonPotentials::LennardJones::Type2>>>;
    
    using LennardJonesType3        = Bond2<addVirial<LennardJones_<CommonPotentials::LennardJones::Type3>>>;
    using LennardJonesType3Const_e = Bond2<addVirial<LennardJonesConst_e_<CommonPotentials::LennardJones::Type3>>>;
    
    using LennardJonesKaranicolasBrooks = Bond2<addVirial<LennardJonesKaranicolasBrooks_>>;
    
    using LennardJonesGaussian = Bond2<addVirial<LennardJonesGaussian_>>;
    using LennardJonesGaussianConst_e_D = Bond2<addVirial<LennardJonesGaussianConst_e_D_>>;

}}}}

#endif
