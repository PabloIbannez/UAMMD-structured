#ifndef __HARMONIC_CONSTRAINT_UNBOUND__
#define __HARMONIC_CONSTRAINT_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{

    struct HarmonicConstraint: public ParameterUpdatable{
        
        std::shared_ptr<ParticleData> pd;

        Box box;

        real K;
        real r0;
            
        struct Parameters{
            
            real K;
            real r0;

        };

        HarmonicConstraint(std::shared_ptr<ParticleData> pd,
                           Parameters par):pd(pd),
                                           K(par.K),
                                           r0(par.r0){}

        ~HarmonicConstraint(){}

        struct forceTransverser{

            real4* force;

            real  K;
            real  r0;

            Box box;

            forceTransverser(real4* force,
                             real   K,
                             real   r0,
                             Box    box):force(force),
                                         K(K),
                                         r0(r0),
                                         box(box){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                    
                if(r2>=r0*r0){
                
                    return make_real4(CommonPotentials::Harmonic::force(rij,r2,K,r0),real(0.0));
                } 

                return make_real4(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     

            return forceTransverser(force,
                                    K,r0,
                                    box);
        }
        
        struct virialTransverser{

            tensor3* virial;

            real  K;
            real  r0;

            Box box;

            virialTransverser(tensor3* virial,
                              real   K,
                              real   r0,
                              Box    box):virial(virial),
                                          K(K),
                                          r0(r0),
                                          box(box){}

            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                return tensor3(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){virial[index_i]+=quantity;}

        };

        virialTransverser getVirialTransverser(){
            
            tensor3* virial = this->pd->getVirial(access::location::gpu, access::mode::readwrite).raw();     

            return virialTransverser(virial,
                                     K,r0,
                                     box);
        }
        
        struct energyTransverser{

            real* energy;

            real  K;
            real  r0;

            Box box;

            energyTransverser(real* energy,
                              real  K,
                              real  r0,
                              Box   box):energy(energy),
                                         K(K),
                                         r0(r0),
                                         box(box){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>r0*r0){

                    return CommonPotentials::Harmonic::energy(rij,r2,K,r0);
                } 

                return real(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy  = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            return energyTransverser(energy,
                                     K,r0,
                                     box);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };

}}}}

#endif
