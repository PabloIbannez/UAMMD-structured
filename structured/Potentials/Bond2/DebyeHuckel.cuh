#ifndef __DEBYE_HUCKEL_BOND2__
#define __DEBYE_HUCKEL_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    template<typename Units_>
    struct DebyeHuckel_{

        struct BondInfo{
            real chgProduct;
            real dielectricConstant;
            real debyeLength;
            real cutOff;
        };
        
        Box box;

        struct Parameters{
            Box box;
        };

        DebyeHuckel_(Parameters par):box(par.box){}

        inline __device__ real3 force(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi){
    
            real3 rij = box.apply_pbc(posj-posi);
            
            real  r2 = dot(rij, rij);

            real3 f=make_real3(0.0);
            if(r2<(bi.cutOff*bi.cutOff)){
                f = CommonPotentials::DebyeHuckel::DebyeHuckel::force<Units_>(rij,r2,
                                                                              bi.chgProduct,
                                                                              bi.dielectricConstant,
                                                                              bi.debyeLength);

                if        (bond_index == i){
                } else if (bond_index == j){
                    f=-f;
                }
            }

            return f;
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const BondInfo &bi) {
            
            real3 rij = box.apply_pbc(posj-posi);
            
            real  r2 = dot(rij, rij);

            real e=real(0.0);
            if(r2<(bi.cutOff*bi.cutOff)){
                e = CommonPotentials::DebyeHuckel::DebyeHuckel::energy<Units_>(rij,r2,
                                                                               bi.chgProduct,
                                                                               bi.dielectricConstant,
                                                                               bi.debyeLength);
            }

            return e/real(2.0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            
            BondInfo bi;
            
            real chgProduct;
            real dielectricConstant;
            real debyeLength;
            real cutOff;
            in>>chgProduct>>dielectricConstant>>debyeLength>>cutOff;
            
            bi.chgProduct = chgProduct;
            bi.dielectricConstant = dielectricConstant;
            bi.debyeLength = debyeLength;
            bi.cutOff = cutOff;
            
            return bi;
        }

    };
    
    template<typename Units_>
    using DebyeHuckel = Bond2<addVirialStress<DebyeHuckel_<Units_>>>;

}}}}

#endif
