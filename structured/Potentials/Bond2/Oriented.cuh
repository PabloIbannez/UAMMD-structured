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
        
        inline __device__ tensor3 virial(int i, int j,
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
    
    struct OrientedHarmonicConst_K_rh_ : public OrientedHarmonic_{
    
        struct BondInfo{};

        real K;
        real rh;
        
        struct Parameters : OrientedHarmonic_::Parameters{
            real K;
            real rh;
        };
        
        OrientedHarmonicConst_K_rh_(Parameters param):OrientedHarmonic_(param),
                                                      K(param.K),
                                                      rh(param.rh){}

        inline __device__ ForceTorque forceTorque(int i, int j,
                                                  int bond_index,
                                                  const real3 &posi,
                                                  const real3 &posj,
                                                  const real4 &diri,
                                                  const real4 &dirj,
                                                  const BondInfo &bi){

            const OrientedHarmonic_::BondInfo biB = {K,rh};

            return OrientedHarmonic_::forceTorque(i,j,bond_index,posi,posj,diri,dirj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const real4 &diri,
                                      const real4 &dirj,
                                      const BondInfo &bi){

            const OrientedHarmonic_::BondInfo biB = {K,rh};

            return OrientedHarmonic_::energy(i,j,bond_index,posi,posj,diri,dirj,biB);
        }
        
        inline __device__ tensor3 virial(int i, int j,
                                         int bond_index,
                                         const real3 &posi,
                                         const real3 &posj,
                                         const real4 &diri,
                                         const real4 &dirj,
                                         const BondInfo &bi) {            
            
            const OrientedHarmonic_::BondInfo biB = {K,rh};
            
            return OrientedHarmonic_::virial(i,j,bond_index,posi,posj,diri,dirj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }
    };
    
    struct OrientedAngular_{

        struct BondInfo{
            real K;
            real theta0;
        };

        Box box;

        struct Parameters{
            Box box;
        };

        OrientedAngular_(Parameters par):box(par.box){}

        inline __device__ ForceTorque forceTorque(int i, int j,
                                                  int bond_index,
                                                  const real3 &posi,
                                                  const real3 &posj,
                                                  const real4 &diri,
                                                  const real4 &dirj,
                                                  const BondInfo &bi){
            
            const real K     = bi.K;
            const real cos0  = cos(bi.theta0);
            
            real3 eix = quaternions::getEx(diri);
            real3 ejx = quaternions::getEx(dirj);
            
            real3 t = -K*(dot(eix,ejx)+cos0)*cross(eix,ejx);

            if        (bond_index == i){
            } else if (bond_index == j){
                t=-t;
            }

            return {make_real4(real(0.0)),
                    make_real4(t,real(0.0))};
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const real4 &diri,
                                      const real4 &dirj,
                                      const BondInfo &bi) {            
            
            const real K     = bi.K;
            const real cos0  = cos(bi.theta0);
            
            real3 eix = quaternions::getEx(diri);
            real3 ejx = quaternions::getEx(dirj);
            
            real costheta = dot(eix,ejx);

            real e = real(0.5)*K*(costheta+cos0)*(costheta+cos0);

            return e/real(2.0);
            
        }
        
        inline __device__ tensor3 virial(int i, int j,
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
            in>>bi.theta0>>bi.K;
            return bi;
        }

    };
    
    struct OrientedAngularConst_K_theta0_ : public OrientedAngular_{
    
        struct BondInfo{};

        real K;
        real theta0;
        
        struct Parameters : OrientedAngular_::Parameters{
            real K;
            real theta0;
        };
        
        OrientedAngularConst_K_theta0_(Parameters param):OrientedAngular_(param),
                                                         K(param.K),
                                                         theta0(param.theta0){}

        inline __device__ ForceTorque forceTorque(int i, int j,
                                                  int bond_index,
                                                  const real3 &posi,
                                                  const real3 &posj,
                                                  const real4 &diri,
                                                  const real4 &dirj,
                                                  const BondInfo &bi){

            const OrientedAngular_::BondInfo biB = {K,theta0};

            return OrientedAngular_::forceTorque(i,j,bond_index,posi,posj,diri,dirj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const real4 &diri,
                                      const real4 &dirj,
                                      const BondInfo &bi){

            const OrientedAngular_::BondInfo biB = {K,theta0};

            return OrientedAngular_::energy(i,j,bond_index,posi,posj,diri,dirj,biB);
        }
        
        inline __device__ tensor3 virial(int i, int j,
                                         int bond_index,
                                         const real3 &posi,
                                         const real3 &posj,
                                         const real4 &diri,
                                         const real4 &dirj,
                                         const BondInfo &bi) {            
            
            const OrientedAngular_::BondInfo biB = {K,theta0};
            
            return OrientedAngular_::virial(i,j,bond_index,posi,posj,diri,dirj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }
    };
    
    struct OrientedDihedral_{

        struct BondInfo{
            real K;
            real phi0;
        };

        Box box;

        struct Parameters{
            Box box;
        };

        OrientedDihedral_(Parameters par):box(par.box){}

        inline __device__ ForceTorque forceTorque(int i, int j,
                                                  int bond_index,
                                                  const real3 &posi,
                                                  const real3 &posj,
                                                  const real4 &diri,
                                                  const real4 &dirj,
                                                  const BondInfo &bi){
            
            const real3 rij = box.apply_pbc(posj-posi);
            
            const real K     = bi.K;
            const real phi0  = bi.phi0;
            
            real3 eiz = quaternions::getEz(diri);
            real3 ejz = quaternions::getEz(dirj);

            real3 F = eiz;
            real3 G = -rij;
            real3 H = ejz;

            real3 A = cross(F,G);
            real3 B = cross(H,G);

            real gn = length(G);
            real an = length(A);
            real bn = length(B);

            real cosphi = dot(A,B)/(an*bn);
            real sinphi = dot(cross(B,A),G)/(an*bn*gn);

            cosphi=std::min(real( 1.0),cosphi);
            cosphi=std::max(real(-1.0),cosphi);
            sinphi=std::min(real( 1.0),sinphi);
            sinphi=std::max(real(-1.0),sinphi);

            real fmod = -K*(sinphi*cos(phi0)-cosphi*sin(phi0));

            real3 f = ((dot(F,G)/(an*an*gn))*A-(dot(H,G)/(bn*bn*gn))*B);

            real3 t = make_real3(0);

            if        (bond_index == i){
                f=-fmod*f;
                t=(fmod*(gn/(an*an)))*cross(eiz,A);
            } else if (bond_index == j){
                f=fmod*f;
                t=-(fmod*(gn/(bn*bn)))*cross(ejz,B);
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
            
            const real K     = bi.K;
            const real phi0  = bi.phi0;
            
            real3 eiz = quaternions::getEz(diri);
            real3 ejz = quaternions::getEz(dirj);

            real3 F = eiz;
            real3 G = -rij;
            real3 H = ejz;

            real3 A = cross(F,G);
            real3 B = cross(H,G);

            real gn = length(G);
            real an = length(A);
            real bn = length(B);

            real cosphi = dot(A,B)/(an*bn);
            real sinphi = dot(cross(B,A),G)/(an*bn*gn);

            cosphi=std::min(real( 1.0),cosphi);
            cosphi=std::max(real(-1.0),cosphi);
            sinphi=std::min(real( 1.0),sinphi);
            sinphi=std::max(real(-1.0),sinphi);
            
            real e = K*(real(1.0)+ cos(phi0)*cosphi + sinphi*sin(phi0));

            return e/real(2.0);
            
        }
        
        inline __device__ tensor3 virial(int i, int j,
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
            in>>bi.phi0>>bi.K;
            return bi;
        }

    };
    
    struct OrientedDihedralConst_K_phi0_ : public OrientedDihedral_{
    
        struct BondInfo{};

        real K;
        real phi0;
        
        struct Parameters : OrientedDihedral_::Parameters{
            real K;
            real phi0;
        };
        
        OrientedDihedralConst_K_phi0_(Parameters param):OrientedDihedral_(param),
                                                        K(param.K),
                                                        phi0(param.phi0){}

        inline __device__ ForceTorque forceTorque(int i, int j,
                                                  int bond_index,
                                                  const real3 &posi,
                                                  const real3 &posj,
                                                  const real4 &diri,
                                                  const real4 &dirj,
                                                  const BondInfo &bi){

            const OrientedDihedral_::BondInfo biB = {K,phi0};

            return OrientedDihedral_::forceTorque(i,j,bond_index,posi,posj,diri,dirj,biB);
        
        }
        
        inline __device__ real energy(int i, int j,
                                      int bond_index,
                                      const real3 &posi,
                                      const real3 &posj,
                                      const real4 &diri,
                                      const real4 &dirj,
                                      const BondInfo &bi){

            const OrientedDihedral_::BondInfo biB = {K,phi0};

            return OrientedDihedral_::energy(i,j,bond_index,posi,posj,diri,dirj,biB);
        }
        
        inline __device__ tensor3 virial(int i, int j,
                                         int bond_index,
                                         const real3 &posi,
                                         const real3 &posj,
                                         const real4 &diri,
                                         const real4 &dirj,
                                         const BondInfo &bi) {            
            
            const OrientedDihedral_::BondInfo biB = {K,phi0};
            
            return OrientedDihedral_::virial(i,j,bond_index,posi,posj,diri,dirj,biB);
        }
        
        static __host__ BondInfo readBond(std::istream &in){
            BondInfo bi;
            return bi;
        }
    };
    
    using OrientedHarmonic           = OrientedBond2<OrientedHarmonic_>;
    using OrientedHarmonicConst_K_rh = OrientedBond2<OrientedHarmonicConst_K_rh_>;

    using OrientedAngular               = OrientedBond2<OrientedAngular_>;
    using OrientedAngularConst_K_theta0 = OrientedBond2<OrientedAngularConst_K_theta0_>;
    
    using OrientedDihedral             = OrientedBond2<OrientedDihedral_>;
    using OrientedDihedralConst_K_phi0 = OrientedBond2<OrientedDihedralConst_K_phi0_>;

}}}}

#endif
