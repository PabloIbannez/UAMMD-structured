#ifndef __DIHEDRAL4_BOND4__
#define __DIHEDRAL4_BOND4__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond4{
            
    struct Dihedral4_{
        
        private:

            Box box;

        public:
            
            struct Parameters{
                Box box;
            };

            Dihedral4_(Parameters par): box(par.box){}

            struct BondInfo{
                real4 K;
                real4 phi0;
            };
            
            struct staticFunctions{

                static inline __device__ real f(const real& cos_dih,const real& sin_dih,const BondInfo& bi){
                    
                    real K[4]  ={bi.K.x,
                                 bi.K.y,
                                 bi.K.z,
                                 bi.K.w};

                    real pha[4]={bi.phi0.x,
                                 bi.phi0.y,
                                 bi.phi0.z,
                                 bi.phi0.w};
                
                    real fmod=0;
                    for(int n=1;n<=4;n++){

                        Dihedral_::BondInfo nbi = {n,K[n-1],pha[n-1]};
                        
                        fmod+=Dihedral_::staticFunctions::f(cos_dih,sin_dih,nbi);
                    }

                    return fmod;
                }
                
                static inline __device__ real e(const real& cos_dih,const real& sin_dih,const BondInfo& bi){
                    
                    real K[4]  ={bi.K.x,
                                 bi.K.y,
                                 bi.K.z,
                                 bi.K.w};

                    real pha[4]={bi.phi0.x,
                                 bi.phi0.y,
                                 bi.phi0.z,
                                 bi.phi0.w};
                
                    real e=0;
                    for(int n=1;n<=4;n++){

                        Dihedral_::BondInfo nbi = {n,K[n-1],pha[n-1]};
                        
                        e+=Dihedral_::staticFunctions::e(cos_dih,sin_dih,nbi);
                    }

                    return e;
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
            
            inline __device__ tensor3 virial(int i, int j, int k, int l,
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
                in >> bi.phi0.x >> bi.phi0.y >> bi.phi0.z >> bi.phi0.w >>
                      bi.K.x    >> bi.K.y    >> bi.K.z    >> bi.K.w    ;
                return bi;
            }

    };
    
    using Dihedral4 = Bond4<Dihedral4_>;

}}}}

#endif
