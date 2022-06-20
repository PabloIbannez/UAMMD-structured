#ifndef __DIHEDRAL_BOND4__
#define __DIHEDRAL_BOND4__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond4{
    
    struct Dihedral_{
        
        private:

            Box box;

        public:
            
            struct Parameters{
                Box box;
            };

            Dihedral_(Parameters par): box(par.box){}

            struct BondInfo{
                int  n;
                real K;
                real phi0;
            };
            
            struct staticFunctions{

                static inline __device__ real f(const real& cos_dih,const real& sin_dih,const BondInfo& bi){
                    
                    const int  n   = bi.n;
                    const real K   = bi.K;
                    const real pha = bi.phi0;
                
                    real cosnt=real(1);
                    real sinnt=real(0);
                    
                    //sinnt cosnt computation
                    real tmp;
                    
                    for(int krot=0;krot<n;krot++){
                        tmp   = cosnt*cos_dih - sinnt*sin_dih;
                        sinnt = sinnt*cos_dih + cosnt*sin_dih;
                        cosnt = tmp;
                    }
                    
                    const real cospha = cos(pha);
                    const real sinpha = sin(pha);
                    
                    return -K*n*(cospha*sinnt-cosnt*sinpha);
                }
                
                static inline __device__ real e(const real& cos_dih,const real& sin_dih,const BondInfo& bi){
                    
                    const int  n   = bi.n;
                    const real K   = bi.K;
                    const real pha = bi.phi0;
                
                    real cosnt=real(1);
                    real sinnt=real(0);
                    
                    real tmp;

                    for(int krot=0;krot<n;krot++){
                        tmp   = cosnt*cos_dih - sinnt*sin_dih;
                        sinnt = sinnt*cos_dih + cosnt*sin_dih;
                        cosnt = tmp;
                    }

                    const real cospha = cos(pha);
                    const real sinpha = sin(pha);

                    return K*(real(1.0)+ cospha*cosnt + sinnt*sinpha);
                }
            };

            inline __device__ real3 force(int i, int j, int k, int l,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const real3 &posl,
                                          const BondInfo &bi)
			{
                real3 fi;
                real3 fjk;
                real3 fl;
                
                dihedralForce<BondInfo,staticFunctions>(posi,posj,posk,posl,
                                                        box,
                                                        bi,
                                                        fi,fjk,fl);
                
				if (bond_index == i){
					return -fi;
                } else if (bond_index == j){
                    return fi-fjk;
                } else if (bond_index == k){
                    return fl+fjk;
                } else if (bond_index == l){
                    return -fl;
                } else{
                }
					
                return make_real3(0);
			}
            
            inline __device__ real virial(int i, int j, int k, int l,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const real3 &posl,
                                          const BondInfo &bi)

			{
				return real(0);
            }
            
            inline __device__ tensor3 stress(int i, int j, int k, int l,
                                             int bond_index,
                                             const real3 &posi,
                                             const real3 &posj,
                                             const real3 &posk,
                                             const real3 &posl,
                                             const BondInfo &bi)

			{
				return tensor3(0);
            }


            inline __device__ real energy(int i, int j, int k, int l,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const real3 &posl,
                                          const BondInfo &bi)

			{
                real e;
                
                dihedralEnergy<BondInfo,staticFunctions>(posi,posj,posk,posl,
                                                         box,
                                                         bi,
                                                         e);
                
                return e/real(4.0);
            }
            
            static BondInfo readBond(std::istream &in){
                BondInfo bi;
                in >> bi.n >> bi.phi0 >> bi.K;
                return bi;
            }

    };
            
    struct DihedralConst_n_K_phi0_ : public Dihedral_{
        
        private:

            int  n;
            real K;
            real phi0;

        public:
            
            struct Parameters : Dihedral_::Parameters{
                int  n;
                real K;
                real phi0;
            };
            
            struct BondInfo{};

            DihedralConst_n_K_phi0_(Parameters par): Dihedral_(par),
                                                     n(par.n),
                                                     K(par.K),
                                                     phi0(par.phi0){}

            inline __device__ real3 force(int i, int j, int k, int l,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const real3 &posl,
                                          const BondInfo &bi)
			{

                const Dihedral_::BondInfo biB = {n,K,phi0};

                return Dihedral_::force(i,j,k,l,bond_index,posi,posj,posk,posl,biB);
			}
            
            inline __device__ real virial(int i, int j, int k, int l,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const real3 &posl,
                                          const BondInfo &bi)

			{
                const Dihedral_::BondInfo biB = {n,K,phi0};

                return Dihedral_::virial(i,j,k,l,bond_index,posi,posj,posk,posl,biB);
            }
            
            inline __device__ tensor3 stress(int i, int j, int k, int l,
                                             int bond_index,
                                             const real3 &posi,
                                             const real3 &posj,
                                             const real3 &posk,
                                             const real3 &posl,
                                             const BondInfo &bi)

			{
                const Dihedral_::BondInfo biB = {n,K,phi0};

                return Dihedral_::stress(i,j,k,l,bond_index,posi,posj,posk,posl,biB);
            }


            inline __device__ real energy(int i, int j, int k, int l,
                                          int bond_index,
                                          const real3 &posi,
                                          const real3 &posj,
                                          const real3 &posk,
                                          const real3 &posl,
                                          const BondInfo &bi)

			{
                const Dihedral_::BondInfo biB = {n,K,phi0};

                return Dihedral_::energy(i,j,k,l,bond_index,posi,posj,posk,posl,biB);
            }
            
            static BondInfo readBond(std::istream &in){
                BondInfo bi;
                return bi;
            }

    };
    
    
    using Dihedral                  = Bond4<Dihedral_>;
    using DihedralConst_n_K_phi0 = Bond4<DihedralConst_n_K_phi0_>;

}}}}

#endif
