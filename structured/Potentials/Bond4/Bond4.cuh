#ifndef __BOND4__
#define __BOND4__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond4{

template <class BondType>
struct Bond4{

    public:
        
        static constexpr int nPart = 4;

        struct Parameters: public BondType::Parameters{};

        struct Bond{

            int i,j,k,l;

            typename BondType::BondInfo bondInfo;

        };
                
        static inline Bond readBond(std::shared_ptr<System> sys,
                                    std::string& line) { 
                
            std::stringstream parser(line);
        
            int i,j,k,l;
            if(!(parser>>i>>j>>k>>l)) {
                sys->log<System::CRITICAL>("Line unreable, %s", line.c_str());
            }
        
            Bond bond;
            
            bond.i = i;
            bond.j = j;
            bond.k = k;
            bond.l = l;
        
            bond.bondInfo = BondType::readBond(parser);
        
            //std::cout << line << std::endl;
        
            return bond;
        
        }

        static inline void registerBond(std::shared_ptr<System>        sys,
                                        std::vector<std::vector<int>>& bondsPerParticle,
                                        int bondIndex,
                                        Bond b){

            //i
            bondsPerParticle[b.i].push_back(bondIndex);
            
            //j
            bondsPerParticle[b.j].push_back(bondIndex);
            
            //k
            bondsPerParticle[b.k].push_back(bondIndex);
            
            //k
            bondsPerParticle[b.l].push_back(bondIndex);
        }
    
    private:

        std::shared_ptr<ParticleData> pd;

        Parameters param;

    public:

        Bond4(std::shared_ptr<ParticleData> pd,
              Parameters param):pd(pd),
                                param(param){}

        struct ForceTransverser: public BondType{

            real4* pos;
            real4* force;

            const int* id2index;

            using resultType=real4;

            ForceTransverser(real4* pos,
                             real4* force,
                             const int* id2index,
                             Parameters param):BondType(param),
                                               pos(pos),
                                               force(force),
                                               id2index(id2index){}

            inline __device__ resultType zero(){return make_real4(0.0);}

            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

            inline __device__ resultType compute(const int index_i,Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];
                const int k = id2index[bond.k];
                const int l = id2index[bond.l];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                real3 posk = make_real3(pos[k]);
                real3 posl = make_real3(pos[l]);

                return make_real4(BondType::force(i,j,k,l,index_i,posi,posj,posk,posl,bond.bondInfo),real(0));

            }

            inline __device__ void set(const int& index_i,resultType& quantity){force[index_i]+=quantity;}

        };

        ForceTransverser getForceTransverser(){

            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real4* force = this->pd->getForce(access::location::gpu, 
                                              access::mode::readwrite).raw();

            const int* id2index = pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(pos,
                                    force,
                                    id2index,
                                    param);
        }
        
        struct VirialTransverser: public BondType{

            real4* pos;
            tensor3* virial;

            const int* id2index;

            using resultType=tensor3;

            VirialTransverser(real4* pos,
                              tensor3* virial,
                              const int* id2index,
                              Parameters param):BondType(param),
                                                pos(pos),
                                                virial(virial),
                                                id2index(id2index){}

            inline __device__ resultType zero(){return tensor3(0.0);}

            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

            inline __device__ resultType compute(const int index_i,Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];
                const int k = id2index[bond.k];
                const int l = id2index[bond.l];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                real3 posk = make_real3(pos[k]);
                real3 posl = make_real3(pos[l]);

                return BondType::virial(i,j,k,l,index_i,posi,posj,posk,posl,bond.bondInfo);

            }

            inline __device__ void set(const int& index_i,resultType& quantity){virial[index_i]+=quantity;}

        };

        VirialTransverser getVirialTransverser(){

            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            tensor3* virial = this->pd->getVirial(access::location::gpu, 
                                              access::mode::readwrite).raw();

            const int* id2index = pd->getIdOrderedIndices(access::location::gpu);

            return VirialTransverser(pos,
                                     virial,
                                     id2index,
                                     param);
        }
        
        struct EnergyTransverser: public BondType{
            
            real4* pos;
            real* energy;

            const int* id2index;

            using resultType=real;

            EnergyTransverser(real4* pos,
                              real* energy,
                              const int* id2index,
                              Parameters param):BondType(param),
                                                pos(pos),
                                                energy(energy),
                                                id2index(id2index){}

            inline __device__ resultType zero(){return real(0.0);}

            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

            inline __device__ resultType compute(const int index_i,Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];
                const int k = id2index[bond.k];
                const int l = id2index[bond.l];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                real3 posk = make_real3(pos[k]);
                real3 posl = make_real3(pos[l]);

                return real(BondType::energy(i,j,k,l,index_i,posi,posj,posk,posl,bond.bondInfo));

            }

            inline __device__ void set(const int& index_i,resultType& quantity){energy[index_i]+=quantity;}

        };

        EnergyTransverser getEnergyTransverser(){

            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real* energy = this->pd->getEnergy(access::location::gpu, 
                                              access::mode::readwrite).raw();

            const int* id2index = pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(pos,
                                     energy,
                                     id2index,
                                     param);
        }
        
        void updateBox(Box box){
            if constexpr (has_box<Parameters>::value) {
                param.box = box;
            }
        }

};

template<typename BondInfo, class staticFunctions>
inline __device__ void dihedralEnergy(const real3 &posi,
                                      const real3 &posj,
                                      const real3 &posk,
                                      const real3 &posl,
                                      const Box   &box,
                                      const BondInfo &bi,
                                      real& e){
				
    const real3 dij = box.apply_pbc(posi - posj);
    const real3 djk = box.apply_pbc(posj - posk);
    const real3 dlk = box.apply_pbc(posl - posk);
    
    const real3 aijk = cross(dij,djk);
    const real3 ajkl = cross(dlk,djk);
    
    const real raijk2=dot(aijk,aijk);
    const real rajkl2=dot(ajkl,ajkl);
    
    const real inv_raijkl = rsqrt(raijk2*rajkl2);
    
    real cos_dih = dot(aijk,ajkl)*inv_raijkl;
    cos_dih=std::min(real( 1.0),cos_dih);
    cos_dih=std::max(real(-1.0),cos_dih);
    
    const real rjk     = sqrt(dot(djk,djk));
    
    real sin_dih = dot(aijk,dlk)*rjk*inv_raijkl;
    sin_dih=std::min(real( 1.0),sin_dih);
    sin_dih=std::max(real(-1.0),sin_dih);

    e = staticFunctions::e(cos_dih,sin_dih,bi);
}

template<typename BondInfo, class staticFunctions>
inline __device__ void dihedralForce(const real3 &posi,
                                     const real3 &posj,
                                     const real3 &posk,
                                     const real3 &posl,
                                     const Box   &box,
                                     const BondInfo &bi,
                                     real3& fi,
                                     real3& fjk,
                                     real3& fl){
				
    const real3 dij = box.apply_pbc(posi - posj);
    const real3 djk = box.apply_pbc(posj - posk);
    const real3 dlk = box.apply_pbc(posl - posk);
    
    const real3 aijk = cross(dij,djk);
    const real3 ajkl = cross(dlk,djk);
    
    const real raijk2=dot(aijk,aijk);
    const real rajkl2=dot(ajkl,ajkl);
    
    const real inv_raijk2=real(1.0)/raijk2;
    const real inv_rajkl2=real(1.0)/rajkl2;
    
    const real inv_raijkl = rsqrt(raijk2*rajkl2);
    
    real cos_dih = dot(aijk,ajkl)*inv_raijkl;
    cos_dih=std::min(real( 1.0),cos_dih);
    cos_dih=std::max(real(-1.0),cos_dih);
    
    const real rjk     = sqrt(dot(djk,djk));
    const real inv_rjk = real(1.0)/rjk;
    
    real sin_dih = dot(aijk,dlk)*rjk*inv_raijkl;
    sin_dih=std::min(real( 1.0),sin_dih);
    sin_dih=std::max(real(-1.0),sin_dih);
   
    //

    const real dot_ijk = dot(dij,djk);
    const real dot_jkl = dot(djk,dlk);
    
    const real3 grad_i  =  rjk*inv_raijk2*aijk;
    const real3 grad_jk = (-dot_ijk*inv_raijk2*aijk+dot_jkl*inv_rajkl2*ajkl)*inv_rjk;
    const real3 grad_l  = -rjk*inv_rajkl2*ajkl;
    
    //

    real fmod = staticFunctions::f(cos_dih,sin_dih,bi);

    fi  = -fmod*grad_i;
    fjk = -fmod*grad_jk;
    fl  = -fmod*grad_l;
}

}}}}

#include "Dihedral.cuh"
#include "Dihedral4.cuh"

#endif
