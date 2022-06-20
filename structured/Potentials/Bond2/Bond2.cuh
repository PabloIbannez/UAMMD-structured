#ifndef __BOND2__
#define __BOND2__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond2{

template <class BondType>
class Bond2Base_: public ParameterUpdatable{

    public:
        
        static constexpr int nPart = 2;

        struct Parameters: public BondType::Parameters{};

        struct Bond{

            int i,j;

            typename BondType::BondInfo bondInfo;

        };
            
        static inline Bond readBond(std::shared_ptr<System> sys,
                                    std::string& line) { 
                
            std::stringstream parser(line);
        
            int i,j;
            if(!(parser>>i>>j)) {
                sys->log<System::CRITICAL>("Line unreable, %s", 
                                            line.c_str());
            }
        
            Bond bond;
            
            bond.i = i;
            bond.j = j;

            bond.bondInfo = BondType::readBond(parser);
        
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
        
        }

    protected:
        
        std::shared_ptr<ParticleData> pd;
        
        Parameters param;

    public:

        Bond2Base_(std::shared_ptr<ParticleData> pd,
                   Parameters param):pd(pd),
                                     param(param){}

        void updateBox(Box box){
            if constexpr (has_box<Parameters>::value) {
                param.box = box;
            }
        }

};

template<class BondType>
class Bond2 : public Bond2Base_<BondType>{

        using Bond2Base = Bond2Base_<BondType>;
        
    public:

        struct Parameters: public Bond2Base::Parameters{};
    
        Bond2(std::shared_ptr<ParticleData> pd,
              Parameters param):Bond2Base(pd,param){}


        struct ForceTransverser: public BondType{

            real4* pos;
            real4* force;
        
            const int* id2index;
            
            using resultType=real4;

            ForceTransverser(real4* pos,
                             real4* force,
                             const int* id2index,
                             typename Bond2Base::Parameters param):BondType(param),
                                                                   pos(pos),
                                                                   force(force),
                                                                   id2index(id2index){}
                
            inline __device__ resultType zero(){return make_real4(0.0);}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);

                return make_real4(BondType::force(i,j,index_i,posi,posj,bond.bondInfo),real(0));
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){force[index_i]+=quantity;}

        };
            
        ForceTransverser getForceTransverser(){
            
            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real4* force = this->pd->getForce(access::location::gpu, 
                                              access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(pos,
                                    force,
                                    id2index,
                                    Bond2Base::param);
        }
        
        struct VirialTransverser: public BondType{

            real4*   pos;
            real* virial;
        
            const int* id2index;
            
            using resultType=real;

            VirialTransverser(real4* pos,
                              real* virial,
                              const int* id2index,
                              typename Bond2Base::Parameters param):BondType(param),
                                                                    pos(pos),
                                                                    virial(virial),
                                                                    id2index(id2index){}
                
            inline __device__ resultType zero(){return real(0.0);}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);

                return BondType::virial(i,j,index_i,posi,posj,bond.bondInfo);
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){virial[index_i]+=quantity;}

        };
            
        VirialTransverser getVirialTransverser(){
            
            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real* virial = this->pd->getVirial(access::location::gpu, 
                                               access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return VirialTransverser(pos,
                                     virial,
                                     id2index,
                                     Bond2Base::param);
        }
        
        struct EnergyTransverser: public BondType{

            real4* pos;
            real* energy;
        
            const int* id2index;
            
            using resultType=real;

            EnergyTransverser(real4* pos,
                              real* energy,
                              const int* id2index,
                              typename Bond2Base::Parameters param):BondType(param),
                                                                    pos(pos),
                                                                    energy(energy),
                                                                    id2index(id2index){}
                
            inline __device__ resultType zero(){return real(0.0);}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);

                return real(BondType::energy(i,j,index_i,posi,posj,bond.bondInfo));
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){energy[index_i]+=quantity;}

        };
            
        EnergyTransverser getEnergyTransverser(){
            
            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real* energy = this->pd->getEnergy(access::location::gpu, 
                                              access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(pos,
                                     energy,
                                     id2index,
                                     Bond2Base::param);
        }
        
        struct StressTransverser: public BondType{

            real4*   pos;
            tensor3* stress;
        
            const int* id2index;
            
            using resultType=tensor3;

            StressTransverser(real4* pos,
                              tensor3* stress,
                              const int* id2index,
                              typename Bond2Base::Parameters param):BondType(param),
                                                                    pos(pos),
                                                                    stress(stress),
                                                                    id2index(id2index){}
                
            inline __device__ resultType zero(){return tensor3(0.0);}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);

                return BondType::stress(i,j,index_i,posi,posj,bond.bondInfo);
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){stress[index_i]+=quantity;}

        };
            
        StressTransverser getStressTransverser(){
            
            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            tensor3* stress = this->pd->getStress(access::location::gpu, 
                                              access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return StressTransverser(pos,
                                     stress,
                                     id2index,
                                     Bond2Base::param);
        }
};

struct ForceTorque{
    real4 force;
    real4 torque;
};

template<class BondType>
class OrientedBond2 : public Bond2Base_<BondType>{

        using Bond2Base = Bond2Base_<BondType>;
        
    public:

        struct Parameters: public Bond2Base::Parameters{};
    
        OrientedBond2(std::shared_ptr<ParticleData> pd,
                      Parameters param):Bond2Base(pd,param){}

        struct ForceTransverser: public BondType{

            real4*      pos;
            real4* dir;
            
            real4* force;
            real4* torque;
        
            const int* id2index;
            
            using resultType=ForceTorque;

            ForceTransverser(real4*      pos,
                             real4* dir,
                             real4* force,
                             real4* torque,
                             const int* id2index,
                             typename Bond2Base::Parameters param):BondType(param),
                                                                   pos(pos),
                                                                   dir(dir),
                                                                   force(force),
                                                                   torque(torque),
                                                                   id2index(id2index){}
                
            inline __device__ resultType zero(){return {.force =make_real4(0.0),
                                                        .torque=make_real4(0.0)};}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total.force +=current.force;
                                                                                          total.torque+=current.torque;}

            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                
                real4 diri = dir[i];
                real4 dirj = dir[j];

                return BondType::forceTorque(i,j,index_i,posi,posj,diri,dirj,bond.bondInfo);
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){force[index_i] +=quantity.force;
                                                                                torque[index_i]+=quantity.torque;}
        };
            
        ForceTransverser getForceTransverser(){
            
            real4*      pos    = this->pd->getPos(access::location::gpu, 
                                             access::mode::read).raw();     
            real4* dir    = this->pd->getDir(access::location::gpu, 
                                             access::mode::read).raw();     
            
            real4* torque = this->pd->getTorque(access::location::gpu, 
                                                access::mode::readwrite).raw();
            real4* force  = this->pd->getForce(access::location::gpu, 
                                               access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(pos,
                                    dir,
                                    force,
                                    torque,
                                    id2index,
                                    Bond2Base::param);
        }
        
        struct VirialTransverser: public BondType{

            real4*   pos;
            real* virial;
            
            real4* dir;
        
            const int* id2index;
            
            using resultType=real;

            VirialTransverser(real4* pos,
                              real4* dir,
                              real* virial,
                              const int* id2index,
                              typename Bond2Base::Parameters param):BondType(param),
                                                                    pos(pos),
                                                                    dir(dir),
                                                                    virial(virial),
                                                                    id2index(id2index){}
                
            inline __device__ resultType zero(){return real(0.0);}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                
                real4 diri = dir[i];
                real4 dirj = dir[j];

                return BondType::virial(i,j,index_i,posi,posj,diri,dirj,bond.bondInfo);
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){virial[index_i]+=quantity;}

        };
            
        VirialTransverser getVirialTransverser(){
            
            real4*      pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real4* dir   = this->pd->getDir(access::location::gpu, 
                                            access::mode::read).raw();     
            
            real* virial = this->pd->getVirial(access::location::gpu, 
                                               access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return VirialTransverser(pos,
                                     dir,
                                     virial,
                                     id2index,
                                     Bond2Base::param);
        }
        
        struct EnergyTransverser: public BondType{

            real4* pos;
            real4* dir;

            real* energy;
        
            const int* id2index;
            
            using resultType=real;

            EnergyTransverser(real4* pos,
                              real4* dir,
                              real* energy,
                              const int* id2index,
                              typename Bond2Base::Parameters param):BondType(param),
                                                                    pos(pos),
                                                                    dir(dir),
                                                                    energy(energy),
                                                                    id2index(id2index){}
                
            inline __device__ resultType zero(){return real(0.0);}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                
                real4 diri = dir[i];
                real4 dirj = dir[j];

                return real(BondType::energy(i,j,index_i,posi,posj,diri,dirj,bond.bondInfo));
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){energy[index_i]+=quantity;}

        };
            
        EnergyTransverser getEnergyTransverser(){
            
            real4* pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real4* dir   = this->pd->getDir(access::location::gpu, 
                                            access::mode::read).raw();     
            real* energy = this->pd->getEnergy(access::location::gpu, 
                                              access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(pos,
                                     dir,
                                     energy,
                                     id2index,
                                     Bond2Base::param);
        }
        
        struct StressTransverser: public BondType{

            real4*   pos;
            tensor3* stress;
            
            real4* dir;
        
            const int* id2index;
            
            using resultType=tensor3;

            StressTransverser(real4* pos,
                              real4* dir,
                              tensor3* stress,
                              const int* id2index,
                              typename Bond2Base::Parameters param):BondType(param),
                                                                    pos(pos),
                                                                    dir(dir),
                                                                    stress(stress),
                                                                    id2index(id2index){}
                
            inline __device__ resultType zero(){return tensor3(0.0);}
                
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,typename Bond2Base::Bond bond){

                const int i = id2index[bond.i];
                const int j = id2index[bond.j];

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                
                real4 diri = dir[i];
                real4 dirj = dir[j];

                return BondType::stress(i,j,index_i,posi,posj,diri,dirj,bond.bondInfo);
                    
            }
                
            inline __device__ void set(const int& index_i,resultType& quantity){stress[index_i]+=quantity;}

        };
            
        StressTransverser getStressTransverser(){
            
            real4*      pos   = this->pd->getPos(access::location::gpu, 
                                            access::mode::read).raw();     
            real4* dir   = this->pd->getDir(access::location::gpu, 
                                            access::mode::read).raw();     
            
            tensor3* stress = this->pd->getStress(access::location::gpu, 
                                              access::mode::readwrite).raw();
            
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return StressTransverser(pos,
                                     dir,
                                     stress,
                                     id2index,
                                     Bond2Base::param);
        }
};


template<class potential>
struct addVirialStress: public potential{

    using potential::potential;

    inline __device__ real virial(int i, int j,
                                  int bond_index,
                                  const real3 &posi,
                                  const real3 &posj,
                                  const typename potential::BondInfo &bi){

        return computeVirial(potential::box.apply_pbc(posj-posi),
                             potential::force(i,j,i,posi,posj,bi)); //Not a typo
    }
    
    inline __device__ tensor3 stress(int i, int j,
                                     int bond_index,
                                     const real3 &posi,
                                     const real3 &posj,
                                     const typename potential::BondInfo &bi){

        return computeStress(potential::box.apply_pbc(posj-posi),
                             potential::force(i,j,i,posi,posj,bi)); //Not a typo
    }

};

}}}}

#include "Harmonic.cuh"
#include "Gaussian.cuh"
#include "Steric.cuh"
#include "Fene.cuh"
#include "Morse.cuh"
#include "LennardJones.cuh"
#include "DebyeHuckel.cuh"

#include "Oriented.cuh"

#endif
