#ifndef __BOND1__
#define __BOND1__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond1{

template <class BondType>
class Bond1{

    public:
        
        static constexpr int nPart = 1;

        struct Parameters: public BondType::Parameters{};

        struct Bond{

            int i;

            typename BondType::BondInfo bondInfo;

        };
            
        static inline Bond readBond(std::shared_ptr<System> sys,
                                    std::string& line) { 
                
            std::stringstream parser(line);
        
            int i;
            if(!(parser>>i)) {
                sys->log<System::CRITICAL>("Line unreable, %s", 
                                            line.c_str());
            }
        
            Bond bond;
            
            bond.i = i;

            bond.bondInfo = BondType::readBond(parser);
        
            return bond;
        
        }

        static inline void registerBond(std::shared_ptr<System>        sys,
                                        std::vector<std::vector<int>>& bondsPerParticle,
                                        int bondIndex,
                                        Bond b){

            //i
            bondsPerParticle[b.i].push_back(bondIndex);
        
        }

    private:
        
        std::shared_ptr<ParticleData> pd;
        
        Parameters param;

    public:

        Bond1(std::shared_ptr<ParticleData> pd,
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

                real3 posi = make_real3(pos[i]);

                return make_real4(BondType::force(i,index_i,posi,bond.bondInfo),real(0));
                    
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

            real4*   pos;
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

                real3 posi = make_real3(pos[i]);

                return BondType::virial(i,index_i,posi,bond.bondInfo);
                    
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

                real3 posi = make_real3(pos[i]);

                return real(BondType::energy(i,index_i,posi,bond.bondInfo));
                    
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

template<class potential>
struct addVirial: public potential{

    using potential::potential;

    inline __device__ tensor3 virial(int i,
                                     int bond_index,
                                     const real3 &posi,
                                     const typename potential::BondInfo &bi){

        return computeVirial(potential::box.apply_pbc(bi.pos-posi),
                             potential::force(i,i,posi,bi)); //Not a typo
    }

};

}}}}

#include "Harmonic.cuh"

#endif
