#ifndef __BOND3__
#define __BOND3__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bond3{

template <class BondType>
class Bond3{

    public:
        
        static constexpr int nPart = 3;

        struct Parameters: public BondType::Parameters{};

        struct Bond{

            int i,j,k;

            typename BondType::BondInfo bondInfo;

        };
                
        static inline Bond readBond(std::shared_ptr<System> sys,
                                    std::string& line) { 
                
            std::stringstream parser(line);
        
            int i,j,k;
            if(!(parser>>i>>j>>k)) {
                sys->log<System::CRITICAL>("Line unreable, %s", line.c_str());
            }
        
            Bond bond;
            
            bond.i = i;
            bond.j = j;
            bond.k = k;
        
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
        }
    
    private:

        std::shared_ptr<ParticleData> pd;

        Parameters param;

    public:

        Bond3(std::shared_ptr<ParticleData> pd,
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

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                real3 posk = make_real3(pos[k]);
    
                return make_real4(BondType::force(i,j,k,index_i,posi,posj,posk,bond.bondInfo),real(0));

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

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                real3 posk = make_real3(pos[k]);

                return BondType::virial(i,j,k,index_i,posi,posj,posk,bond.bondInfo);

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

                real3 posi = make_real3(pos[i]);
                real3 posj = make_real3(pos[j]);
                real3 posk = make_real3(pos[k]);
    
                return real(BondType::energy(i,j,k,index_i,posi,posj,posk,bond.bondInfo));

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
inline __device__ void angularEnergy(const real3 &posi,
                                     const real3 &posj,
                                     const real3 &posk,
                                     const Box   &box,
                                     const BondInfo &bi,
                                     real& e){
                
    //         i -------- j -------- k
    //             <- rji     rjk ->
    //Compute distances and vectors
    //---rji---
    const real3 rji =  box.apply_pbc(posi - posj);
    const real rji2 = dot(rji, rji);
    //---rkj---
    const real3 rjk =  box.apply_pbc(posk - posj);
    const real rjk2 = dot(rjk, rjk);

    const real inv_rjirjk = rsqrt(rji2*rjk2);
    const real inv_rji2 = real(1.0)/rji2;
    const real inv_rjk2 = real(1.0)/rjk2;
    
    real cijk = dot(rji, rjk)*inv_rjirjk; //cijk = cos (theta) = rji*rkj / mod(rji)*mod(rkj)
    //Cos must stay in range
    cijk = min(real( 1.0),cijk);
    cijk = max(real(-1.0),cijk);

    const real ang = acos(cijk);

    e = staticFunctions::e(ang,bi);     
}

template<typename BondInfo, class staticFunctions>
inline __device__ void angularForce(const real3 &posi,
                                    const real3 &posj,
                                    const real3 &posk,
                                    const Box   &box,
                                    const BondInfo &bi,
                                    real3& fi,
                                    real3& fk){
                
    //         i -------- j -------- k
    //             <- rji     rjk ->
    //Compute distances and vectors
    //---rji---
    const real3 rji =  box.apply_pbc(posi - posj);
    const real rji2 = dot(rji, rji);
    //---rkj---
    const real3 rjk =  box.apply_pbc(posk - posj);
    const real rjk2 = dot(rjk, rjk);

    const real inv_rjirjk = rsqrt(rji2*rjk2);
    const real inv_rji2 = real(1.0)/rji2;
    const real inv_rjk2 = real(1.0)/rjk2;
    
    real cijk = dot(rji, rjk)*inv_rjirjk; //cijk = cos (theta) = rji*rkj / mod(rji)*mod(rkj)
    //Cos must stay in range
    cijk = min(real( 1.0),cijk);
    cijk = max(real(-1.0),cijk);

    const real ang = acos(cijk);

    real fmod = staticFunctions::f(ang,bi); 

    real sijk = sin(ang);
    //Sin must be bigger than 0
    //sijk = max(std::numeric_limits<real>::min(),sijk);
    sijk = max(real(1e-6),sijk);
    fmod=fmod/sijk;

    const real crji = cijk*inv_rji2;
    const real crjk = cijk*inv_rjk2;

    fi = fmod*(rjk*inv_rjirjk-rji*crji);
    fk = fmod*(rji*inv_rjirjk-rjk*crjk);
    
    /*printf("i:%i j:%i k:%i \\
            rji_x:%f rji_y:%f rji_z:%f \\ 
            rjk_x:%f rjk_y:%f rjk_z:%f \\
            angle:%f angle0%f anglediff:%f \\ 
            K:%f \\
            fmod:%f \\
            fi%i:%f %f %f \\ 
            fj%i:%f %f %f\n",
            i,j,k,
            rji.x,rji.y,rji.z,
            rjk.x,rjk.y,rjk.z,
            ang,bi.ang0,ang-bi.ang0,
            bi.K,
            fmod,
            i,-fi.x,-fi.y,-fi.z,
            k,-fk.x,-fk.y,-fk.z);*/
               
}


}}}}

#include "HarmonicAngular.cuh"
#include "KratkyPorod.cuh"
#include "BestChenHummerAngular.cuh"

#endif
