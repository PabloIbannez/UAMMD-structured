#ifndef __STERIC_UNBOUND__
#define __STERIC_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{

    template<int power>
    struct StericPBC: public ParameterUpdatable{
        
        std::shared_ptr<ParticleData> pd;

        Box box;

        real stericCutOff;
            
        bool useInnerRadius;
            
        struct Parameters{

            real stericCutOff;

            bool useInnerRadius = false;
        };

        StericPBC(std::shared_ptr<ParticleData> pd,
                  Parameters par):pd(pd),
                                  stericCutOff(par.stericCutOff),
                                  useInnerRadius(par.useInnerRadius){}

        ~StericPBC(){}

        struct forceTransverser{

            real4* force;

            real* radius;
            real* epsilon;

            Box box;

            real stericCutOff2;

            forceTransverser(real4* force,
                             real*  radius,
                             real*  epsilon,
                             Box    box,
                             real   stericCutOff2):force(force),
                                                   radius(radius),
                                                   epsilon(epsilon),
                                                   box(box),
                                                   stericCutOff2(stericCutOff2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=stericCutOff2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = sqrt(epsilon[index_i]*epsilon[index_j]);

                    return make_real4(CommonPotentials::Steric::Steric::force<power>(rij,r2,eps,sigma),real(0.0));
                } 

                return make_real4(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius;
            if(useInnerRadius){
                radius = this->pd->getInnerRadius(access::location::gpu, access::mode::readwrite).raw();     
            } else {
                radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            }

            real* epsilon = this->pd->getEpsilon(access::location::gpu, access::mode::readwrite).raw();     

            return forceTransverser(force,
                                    radius,
                                    epsilon,
                                    box,
                                    stericCutOff*stericCutOff);
        }
        
        struct virialTransverser{

            tensor3* virial;

            real* radius;
            real* epsilon;

            Box box;

            real stericCutOff2;

            virialTransverser(tensor3* virial,
                              real*   radius,
                              real*   epsilon,
                              Box     box,
                              real    stericCutOff2):virial(virial),
                                                     radius(radius),
                                                     epsilon(epsilon),
                                                     box(box),
                                                     stericCutOff2(stericCutOff2){}

            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=stericCutOff2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = sqrt(epsilon[index_i]*epsilon[index_j]);
                    
                    return CommonPotentials::Steric::Steric::virial<power>(rij,r2,eps,sigma);
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

            real* epsilon = this->pd->getEpsilon(access::location::gpu, access::mode::readwrite).raw();     

            return virialTransverser(virial,
                                     radius,
                                     epsilon,
                                     box,
                                     stericCutOff*stericCutOff);
        }
        
        struct energyTransverser{

            real* energy;

            real* radius;
            real* epsilon;

            Box box;

            real stericCutOff2;

            energyTransverser(real* energy,
                              real*   radius,
                              real*   epsilon,
                              Box     box,
                              real    stericCutOff2):energy(energy),
                                                     radius(radius),
                                                     epsilon(epsilon),
                                                     box(box),
                                                     stericCutOff2(stericCutOff2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=stericCutOff2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = sqrt(epsilon[index_i]*epsilon[index_j]);

                    return CommonPotentials::Steric::Steric::energy<power>(rij,r2,eps,sigma);
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

            real* epsilon = this->pd->getEpsilon(access::location::gpu, access::mode::readwrite).raw();     

            return energyTransverser(energy,
                                     radius,
                                     epsilon,
                                     box,
                                     stericCutOff*stericCutOff);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }
    };


    template<int power>
    struct StericConstantPBC: public ParameterUpdatable{
        
        std::shared_ptr<ParticleData> pd;

        Box box;

        real stericCutOff;
            
        real sigma;
        real epsilon;
            
        struct Parameters{
            
            real stericCutOff;

            real sigma;
            real epsilon;
        };

        StericConstantPBC(std::shared_ptr<ParticleData> pd,
                          Parameters par):pd(pd),
                                          stericCutOff(par.stericCutOff),
                                          sigma(par.sigma),
                                          epsilon(par.epsilon){}

        ~StericConstantPBC(){}

        struct forceTransverser{

            real4* force;

            real  sigma;
            real  epsilon;

            Box box;

            real stericCutOff2;

            forceTransverser(real4* force,
                             real   sigma,
                             real   epsilon,
                             Box    box,
                             real   stericCutOff2):force(force),
                                                   sigma(sigma),
                                                   epsilon(epsilon),
                                                   box(box),
                                                   stericCutOff2(stericCutOff2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                    
                if(r2>0 and r2<=stericCutOff2){

                    return make_real4(CommonPotentials::Steric::Steric::force<power>(rij,r2,epsilon,sigma),real(0.0));
                } 

                return make_real4(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     

            return forceTransverser(force,
                                    sigma,
                                    epsilon,
                                    box,
                                    stericCutOff*stericCutOff);
        }
        
        struct virialTransverser{

            tensor3* virial;

            real  sigma;
            real  epsilon;

            Box box;

            real stericCutOff2;

            virialTransverser(tensor3* virial,
                              real   sigma,
                              real   epsilon,
                              Box    box,
                              real   stericCutOff2):virial(virial),
                                                    sigma(sigma),
                                                    epsilon(epsilon),
                                                    box(box),
                                                    stericCutOff2(stericCutOff2){}

            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=stericCutOff2){

                    return CommonPotentials::Steric::Steric::virial<power>(rij,r2,epsilon,sigma);
                } 

                return tensor3(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){virial[index_i]+=quantity;}

        };

        virialTransverser getVirialTransverser(){
            
            tensor3* virial = this->pd->getVirial(access::location::gpu, access::mode::readwrite).raw();     

            return virialTransverser(virial,
                                     sigma,
                                     epsilon,
                                     box,
                                     stericCutOff*stericCutOff);
        }
        
        struct energyTransverser{

            real* energy;

            real  sigma;
            real  epsilon;

            Box box;

            real stericCutOff2;

            energyTransverser(real* energy,
                              real  sigma,
                              real  epsilon,
                              Box   box,
                              real  stericCutOff2):energy(energy),
                                                   sigma(sigma),
                                                   epsilon(epsilon),
                                                   box(box),
                                                   stericCutOff2(stericCutOff2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=stericCutOff2){

                    return CommonPotentials::Steric::Steric::energy<power>(rij,r2,epsilon,sigma);
                } 

                return real(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy  = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            return energyTransverser(energy,
                                     sigma,
                                     epsilon,
                                     box,
                                     stericCutOff*stericCutOff);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };

}}}}

#endif
