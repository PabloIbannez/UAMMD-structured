#ifndef __ORIENTED_BOND2__
#define __ORIENTED_BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

    struct OrientedHarmonic_{

        struct BondInfo{
            real K;
            real rh;
        };

        Box box;

        struct Parameters{
            Box box;
        };

        OrientedHarmonic_(Parameters par):box(par.box){}

        inline __device__ ForceTorque forceTorque(int i, int j,
                                                  int bond_index,
                                                  const real3 &posi,
                                                  const real3 &posj,
                                                  const real4 &diri,
                                                  const real4 &dirj,
                                                  const BondInfo &bi){
            
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real K   = bi.K;
            const real rh  = bi.rh;
            
            real3 eix = quaternions::getEx(diri);
            real3 ejx = quaternions::getEx(dirj);
            
            real3 dr = rij-rh*(eix+ejx);

            real3 f = K*dr;
            real3 t = make_real3(0.0);

            if        (bond_index == i){
                t=K*rh*cross(eix,dr);
            } else if (bond_index == j){
                t=K*rh*cross(ejx,dr);
                f=-f;
            }

            return {make_real4(f,real(0.0)),
                    make_real4(t,real(0.0))};
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const real4 &diri,
                                      const real4 &dirj,
                                      const BondInfo &bi) {            
            
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real K   = bi.K;
            const real rh  = bi.rh;
            
            real3 eix = quaternions::getEx(diri);
            real3 ejx = quaternions::getEx(dirj);
            
            real3 dr = rij-rh*(eix+ejx);

            real e = real(0.5)*K*dot(dr,dr);

            return e/real(2.0);
        }
        
        inline __device__ real virial(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const real4 &diri,
                                      const real4 &dirj,
                                      const BondInfo &bi) {            
            
            return real(0);
        }
        
        inline __device__ tensor3 stress(int i, int j,
                                         int bond_index,
                                         const real3 &posi,
                                         const real3 &posj,
                                         const real4 &diri,
                                         const real4 &dirj,
                                         const BondInfo &bi) {            
            
            return tensor3(0);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            in>>bi.rh>>bi.K;
            return bi;
        }

    };
    
    using OrientedHarmonic           = OrientedBond2<OrientedHarmonic_>;

}}}}

#endif
