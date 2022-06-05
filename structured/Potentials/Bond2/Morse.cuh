#ifndef __MORSE_BOND2__
#define __MORSE_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    struct Morse_{

        struct BondInfo{
            real r0;
            real E;
            real D;
        };
        
        Box box;

        struct Parameters{
            Box box;
        };
        
        Morse_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const real3 rij = box.apply_pbc(posj-posi);
            
            const real r0 = bi.r0;
            const real E  = bi.E;
            const real D  = bi.D;
            
            const real r2 = dot(rij, rij);

            real3 f = CommonPotentials::Morse::force(rij,r2,E,r0,D);
            
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
                                      const BondInfo &bi){
            
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real r0 = bi.r0;
            const real E  = bi.E;
            const real D  = bi.D;
            
            const real r2 = dot(rij, rij);

            real e = CommonPotentials::Morse::energy(rij,r2,E,r0,D);
            
            return e;
        }
        
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.r0 >> bi.E >> bi.D;
            return bi;
        }
    };
    
    struct MorseConst_D_ : public Morse_{

        struct BondInfo{
            real r0;
            real E;
        };
        
        struct Parameters:public Morse_::Parameters{
            real D;
        };
            
        real D;
        
        MorseConst_D_(Parameters par):Morse_(par),
                                      D(par.D){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Morse_::BondInfo biB = {.r0=bi.r0,.E=bi.E,.D=D};

            return Morse_::force(i,j,bond_index,posi,posj,biB);
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
            
            const Morse_::BondInfo biB = {.r0=bi.r0,.E=bi.E,.D=D};

            return Morse_::energy(i,j,bond_index,posi,posj,biB);
        }
        
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.r0 >> bi.E;
            return bi;
        }
    };
    
    struct MorseConst_r0_E_D_ : public Morse_{

        struct BondInfo{};
        
        struct Parameters:public Morse_::Parameters{
            real r0;
            real E;
            real D;
        };
            
        real r0;
        real E;
        real D;
        
        MorseConst_r0_E_D_(Parameters par):Morse_(par),
                                           r0(par.r0),
                                           E(par.E),
                                           D(par.D){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const Morse_::BondInfo biB = {.r0=r0,.E=E,.D=D};

            return Morse_::force(i,j,bond_index,posi,posj,biB);
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
            
            const Morse_::BondInfo biB = {.r0=r0,.E=E,.D=D};

            return Morse_::energy(i,j,bond_index,posi,posj,biB);
        }
        
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }
    };

    struct MorseWCA_{

        struct BondInfo{
            real r0;
            real E;
            real D;
        };
        
        real eps0;
        
        Box box;

        struct Parameters{
            real eps0;
            
            Box box;
        };
        
        MorseWCA_(Parameters par):eps0(par.eps0),
                                  box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){

            const real3 rij = box.apply_pbc(posj-posi);
            
            const real r0 = bi.r0;
            const real E  = bi.E;
            const real D  = bi.D;
            
            const real r2 = dot(rij, rij);

            real3 f = CommonPotentials::Morse::force(rij,r2,E,r0,D) +
                      CommonPotentials::WCA::Type2::force(rij,r2,eps0,r0);
            
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
                                      const BondInfo &bi){
            
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real r0 = bi.r0;
            const real E  = bi.E;
            const real D  = bi.D;
            
            const real r2 = dot(rij, rij);

            real e = CommonPotentials::Morse::energy(rij,r2,E,r0,D) + 
                     CommonPotentials::WCA::Type2::energy(rij,r2,eps0,r0);
            
            return e;
        }
        
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.r0 >> bi.E >> bi.D;
            return bi;
        }
    };
    
    
    using Morse             = Bond2<addVirial<Morse_>>;
    using MorseConst_D      = Bond2<addVirial<MorseConst_D_>>;
    using MorseConst_r0_E_D = Bond2<addVirial<MorseConst_r0_E_D_>>;
    
    using MorseWCA          = Bond2<addVirial<MorseWCA_>>;

}}}}

#endif
