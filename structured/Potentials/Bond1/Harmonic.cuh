#ifndef __HARMONIC_BOND1__
#define __HARMONIC_BOND1__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond1{

    struct HarmonicAnisotropic_{

        struct BondInfo{
            real3 K;
            real3 r0;
            real3 pos;
        };
        
        Box box;

        struct Parameters{
            Box box;
        };
        
        HarmonicAnisotropic_(Parameters param):box(param.box){}

        inline __device__ real3 force(int i,
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){
            
            const real3 posj = bi.pos;
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real3 K   = bi.K;
            const real3 r0  = bi.r0;

            return CommonPotentials::HarmonicAnisotropic::force(rij,K,r0);
        }
        
        inline __device__ real energy(int i,
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){
            const real3 posj = bi.pos;
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real3 K   = bi.K;
            const real3 r0  = bi.r0;

            const real e = CommonPotentials::HarmonicAnisotropic::energy(rij,K,r0);
            
            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.K.x   >> bi.K.y   >> bi.K.z   
               >> bi.r0.x  >> bi.r0.y  >> bi.r0.z 
               >> bi.pos.x >> bi.pos.y >> bi.pos.z;
            return bi;
        }
    };
    
    struct HarmonicAnisotropicConst_r0_ : public HarmonicAnisotropic_{
    
        struct BondInfo{
            real3 pos;
            real3  K;
        };
            
        real3  r0;

        struct Parameters : HarmonicAnisotropic_::Parameters{
            real3  r0;
        };
        
        HarmonicAnisotropicConst_r0_(Parameters param):HarmonicAnisotropic_(param),
                                                       r0(param.r0){}

        inline __device__ real3 force(int i,                                      
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = { bi.K,
                                                        {r0.x,r0.y,r0.z},
                                                         bi.pos};

            return HarmonicAnisotropic_::force(i,bond_index,posi,biB);
        
        }
        
        inline __device__ real energy(int i,
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = { bi.K,
                                                        {r0.x,r0.y,r0.z},
                                                         bi.pos};

            return HarmonicAnisotropic_::energy(i,bond_index,posi,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.K.x   >> bi.K.y   >> bi.K.z   
               >> bi.pos.x >> bi.pos.y >> bi.pos.z;
            return bi;
        }
    };
    
    struct HarmonicAnisotropicConst_K_r0_ : public HarmonicAnisotropic_{
    
        struct BondInfo{
            real3 pos;
        };
            
        real3  K;
        real3  r0;

        struct Parameters : HarmonicAnisotropic_::Parameters{
            real3  K;
            real3  r0;
        };
        
        HarmonicAnisotropicConst_K_r0_(Parameters param):HarmonicAnisotropic_(param),
                                                         K(param.K),r0(param.r0){}

        inline __device__ real3 force(int i,                                      
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = {{K.x ,K.y  ,K.z},
                                                        {r0.x,r0.y,r0.z},
                                                         bi.pos};

            return HarmonicAnisotropic_::force(i,bond_index,posi,biB);
        
        }
        
        inline __device__ real energy(int i,
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = {{K.x ,K.y  ,K.z},
                                                        {r0.x,r0.y,r0.z},
                                                         bi.pos};

            return HarmonicAnisotropic_::energy(i,bond_index,posi,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.pos.x >> bi.pos.y >> bi.pos.z;
            return bi;
        }
    };
    
    struct HarmonicAnisotropicCommon_K_r0_ : public HarmonicAnisotropic_{
    
        struct BondInfo{
            real  K;
            real  r0;
            real3 pos;
        };

        struct Parameters : HarmonicAnisotropic_::Parameters{};
        
        HarmonicAnisotropicCommon_K_r0_(Parameters param):HarmonicAnisotropic_(param){}

        inline __device__ real3 force(int i,                                      
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = {{bi.K,bi.K,bi.K},
                                                        {bi.r0,bi.r0,bi.r0},
                                                         bi.pos};

            return HarmonicAnisotropic_::force(i,bond_index,posi,biB);
        
        }
        
        inline __device__ real energy(int i,
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = {{bi.K,bi.K,bi.K},
                                                        {bi.r0,bi.r0,bi.r0},
                                                         bi.pos};

            return HarmonicAnisotropic_::energy(i,bond_index,posi,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.K >> bi.r0 
               >> bi.pos.x >> bi.pos.y >> bi.pos.z;
            return bi;
        }
    };
    
    struct HarmonicAnisotropicConstCommon_K_r0_ : public HarmonicAnisotropic_{
    
        struct BondInfo{
            real3 pos;
        };

        struct Parameters : HarmonicAnisotropic_::Parameters{
            real  K;
            real  r0;
        };
            
        real  K;
        real  r0;
        
        HarmonicAnisotropicConstCommon_K_r0_(Parameters param):HarmonicAnisotropic_(param),
                                                               K(param.K),
                                                               r0(param.r0){}

        inline __device__ real3 force(int i,                                      
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = {{K,K,K},
                                                        {r0,r0,r0},
                                                         bi.pos};

            return HarmonicAnisotropic_::force(i,bond_index,posi,biB);
        
        }
        
        inline __device__ real energy(int i,
                                      int bond_index,
                                      const real3 &posi,
                                      const BondInfo &bi){

            const HarmonicAnisotropic_::BondInfo biB = {{K,K,K},
                                                        {r0,r0,r0},
                                                         bi.pos};

            return HarmonicAnisotropic_::energy(i,bond_index,posi,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in >> bi.pos.x >> bi.pos.y >> bi.pos.z;
            return bi;
        }
    };

    using Harmonic                 = Bond1<addVirial<HarmonicAnisotropic_>>;
    using HarmonicConst_r0         = Bond1<addVirial<HarmonicAnisotropicConst_r0_>>;
    using HarmonicConst_K_r0       = Bond1<addVirial<HarmonicAnisotropicConst_K_r0_>>;
    using HarmonicCommon_K_r0      = Bond1<addVirial<HarmonicAnisotropicCommon_K_r0_>>;
    using HarmonicConstCommon_K_r0 = Bond1<addVirial<HarmonicAnisotropicConstCommon_K_r0_>>;

}}}}

#endif
