#ifndef __WCA_POTENTIAL__
#define __WCA_POTENTIAL__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{

    struct WCA: public ParameterUpdatable{
        
        std::shared_ptr<ParticleData> pd;

        Box box;
            
        real epsilon;

        real cutOff;
        
        bool useInnerRadius;
            
        struct Parameters{

            real epsilon;

            real cutOff;
            
            bool useInnerRadius = false;
        };

        WCA(std::shared_ptr<ParticleData> pd,
            Parameters par):pd(pd),
                            epsilon(par.epsilon),
                            cutOff(par.cutOff),
                            useInnerRadius(par.useInnerRadius){}


        ~WCA(){}

        struct forceTransverser{

            real4* force;
            
            real* radius;

            real epsilon;

            Box box;

            real cutOff2;

            forceTransverser(real4* force,
                             real* radius,
                             real epsilon,
                             Box box,
                             real cutOff2):force(force),
                                           radius(radius),
                                           epsilon(epsilon),
                                           box(box),
                                           cutOff2(cutOff2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                if(r2>0 and r2<=cutOff2){

                    const real sigma      = radius[index_i]+radius[index_j];

                    return make_real4(CommonPotentials::WCA::Type2::force(rij,r2,epsilon,sigma),real(0.0));
                } 

                return make_real4(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            real* radius;
            if(useInnerRadius){
                radius = this->pd->getInnerRadius(access::location::gpu, access::mode::readwrite).raw();     
            } else {
                radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            }

            return forceTransverser(force,
                                    radius,
                                    epsilon,
                                    box,
                                    cutOff*cutOff);
        }
        
        struct virialTransverser{

            tensor3* virial;
                              
            real* radius;

            real epsilon;

            Box box;

            real cutOff2;

            virialTransverser(tensor3* virial,
                              real* radius,
                              real epsilon,
                              Box box,
                              real cutOff2):virial(virial),
                                            radius(radius),
                                            epsilon(epsilon),
                                            box(box),
                                            cutOff2(cutOff2){}

            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                    
                if(r2>0 and r2<=cutOff2){
                    
                    const real sigma      = radius[index_i]+radius[index_j];

                    return (CommonPotentials::WCA::Type2::virial(rij,r2,epsilon,sigma));
                } 

                return tensor3(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){virial[index_i]+=quantity;}

        };

        virialTransverser getVirialTransverser(){
            
            tensor3* virial = this->pd->getVirial(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius;
            if(useInnerRadius){
                radius = this->pd->getInnerRadius(access::location::gpu, access::mode::readwrite).raw();     
            } else {
                radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            }

            return virialTransverser(virial,
                                     radius,
                                     epsilon,
                                     box,
                                     cutOff*cutOff);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* radius;

            real epsilon;

            Box box;

            real cutOff2;

            energyTransverser(real* energy,
                              real* radius,
                              real epsilon,
                              Box box,
                              real cutOff2):energy(energy),
                                            radius(radius),
                                            epsilon(epsilon),
                                            box(box),
                                            cutOff2(cutOff2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=cutOff2){
                    
                    const real sigma      = radius[index_i]+radius[index_j];

                    return (CommonPotentials::WCA::Type2::energy(rij,r2,epsilon,sigma));
                } 

                return real(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy  = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius;
            if(useInnerRadius){
                radius = this->pd->getInnerRadius(access::location::gpu, access::mode::readwrite).raw();     
            } else {
                radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            }
            
            return energyTransverser(energy,
                                     radius,
                                     epsilon,
                                     box,
                                     cutOff*cutOff);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };

}}}}

#endif
