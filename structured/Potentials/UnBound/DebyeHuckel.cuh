#ifndef __DH_UNBOUND__
#define __DH_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{

    template<class Topology>
    struct DebyeHuckel_: public ParameterUpdatable{
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<Topology>     top;

        Box box;
            
        real dielectricConstant; 
        real debyeLength;

        real cutOffDst;
            
        struct Parameters{

            real dielectricConstant; 
            real debyeLength;

            real cutOffDst;
        };

        DebyeHuckel_(std::shared_ptr<ParticleGroup> pg,
                     std::shared_ptr<Topology>     top,
                     Parameters par):pg(pg),
                                     pd(pg->getParticleData()),
                                     sys(pg->getParticleData()->getSystem()),
                                     dielectricConstant(par.dielectricConstant),
                                     debyeLength(par.debyeLength),
                                     cutOffDst(par.cutOffDst){}
        
        real getCutOffDst(){
            return cutOffDst;
        }

        ~DebyeHuckel_(){}

        struct forceTransverser{

            real4* force;
            
            real* charge;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOffDst2;

            forceTransverser(real4* force,
                             real* charge,
                             real dielectricConstant, 
                             real debyeLength,
                             Box box,
                             real cutOffDst2):force(force),
                                           charge(charge),
                                           dielectricConstant(dielectricConstant),
                                           debyeLength(debyeLength),
                                           box(box),
                                           cutOffDst2(cutOffDst2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                if(r2>0 and r2<=cutOffDst2){

                    const real chgProduct = charge[index_i]*charge[index_j];

                    return make_real4(CommonPotentials::DebyeHuckel::DebyeHuckel::force<Topology::Units>(rij,r2,chgProduct,
                                                                                                         dielectricConstant,debyeLength),real(0.0));
                } 

                return make_real4(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     

            return forceTransverser(force,
                                    charge,
                                    dielectricConstant, 
                                    debyeLength,
                                    box,
                                    cutOffDst*cutOffDst);
        }
        
        
        struct energyTransverser{

            real* energy;
            
            real* charge;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOffDst2;

            energyTransverser(real* energy,
                              real* charge,
                              real dielectricConstant, 
                              real debyeLength,
                              Box box,
                              real cutOffDst2):energy(energy),
                                            charge(charge),
                                            dielectricConstant(dielectricConstant),
                                            debyeLength(debyeLength),
                                            box(box),
                                            cutOffDst2(cutOffDst2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=cutOffDst2){
                    
                    const real chgProduct = charge[index_i]*charge[index_j];

                    return (CommonPotentials::DebyeHuckel::DebyeHuckel::energy<Topology::Units>(rij,r2,chgProduct,
                                                                                                dielectricConstant,debyeLength));
                } 

                return real(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy  = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
            
            return energyTransverser(energy,
                                     charge,
                                     dielectricConstant, 
                                     debyeLength,
                                     box,
                                     cutOffDst*cutOffDst);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };
    
    template<class Topology>
    using DebyeHuckel = DebyeHuckel_<Topology>;
    
    template<class Topology>
    struct DebyeHuckelSpheres_: public ParameterUpdatable{
        
        std::shared_ptr<System>        sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<Topology>     top;

        Box box;
            
        real dielectricConstant; 
        real debyeLength;

        real cutOffDst;
            
        struct Parameters{

            real dielectricConstant; 
            real debyeLength;

            real cutOffDst;
        };

        DebyeHuckelSpheres_(std::shared_ptr<ParticleGroup> pg,
                            std::shared_ptr<Topology>     top,
                            Parameters par):pg(pg),
                                            pd(pg->getParticleData()),
                                            sys(pg->getParticleData()->getSystem()),
                                            dielectricConstant(par.dielectricConstant),
                                            debyeLength(par.debyeLength),
                                            cutOffDst(par.cutOffDst){}
        
        real getCutOffDst(){
            return cutOffDst;
        }

        ~DebyeHuckelSpheres_(){}

        struct forceTransverser{

            real4* force;
            
            real* radius;
            real* charge;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOffDst2;

            forceTransverser(real4* force,
                             real* radius,
                             real* charge,
                             real dielectricConstant, 
                             real debyeLength,
                             Box box,
                             real cutOffDst2):force(force),
                                           radius(radius),
                                           charge(charge),
                                           dielectricConstant(dielectricConstant),
                                           debyeLength(debyeLength),
                                           box(box),
                                           cutOffDst2(cutOffDst2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                if(r2>0 and r2<=cutOffDst2){

                    const real chgProduct = charge[index_i]*charge[index_j];

                    return make_real4(CommonPotentials::DebyeHuckel::DebyeHuckelSpheres::force<Topology::Units>(rij,r2,chgProduct,
                                                                                                                radius[index_i],radius[index_j],
                                                                                                                dielectricConstant,debyeLength),real(0.0));
                } 

                return make_real4(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            real*  radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     

            return forceTransverser(force,
                                    radius,
                                    charge,
                                    dielectricConstant, 
                                    debyeLength,
                                    box,
                                    cutOffDst*cutOffDst);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* radius;
            real* charge;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOffDst2;

            energyTransverser(real* energy,
                              real* radius,
                              real* charge,
                              real dielectricConstant, 
                              real debyeLength,
                              Box box,
                              real cutOffDst2):energy(energy),
                                            radius(radius),
                                            charge(charge),
                                            dielectricConstant(dielectricConstant),
                                            debyeLength(debyeLength),
                                            box(box),
                                            cutOffDst2(cutOffDst2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=cutOffDst2){
                    
                    const real chgProduct = charge[index_i]*charge[index_j];

                    return (CommonPotentials::DebyeHuckel::DebyeHuckelSpheres::energy<Topology::Units>(rij,r2,chgProduct,
                                                                                                       radius[index_i],radius[index_j],
                                                                                                       dielectricConstant,debyeLength));
                } 

                return real(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy  = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            real*  radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
            
            return energyTransverser(energy,
                                     radius,
                                     charge,
                                     dielectricConstant, 
                                     debyeLength,
                                     box,
                                     cutOffDst*cutOffDst);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };
    
    template<class Topology>
    using DebyeHuckelSpheres = DebyeHuckelSpheres_<Topology>;
    
    template<class Topology>
    struct DebyeHuckelDistanceDependentDielectric_: public ParameterUpdatable{
            
        static constexpr real wd  = real(78.4);
        static constexpr real A   = real(-8.5525);
        static constexpr real B   = wd-A;
        static constexpr real k   = real(7.7839);
        static constexpr real lmd = real(0.003627);
            
        std::shared_ptr<System>        sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<Topology>     top;

        Box box;
            
        real debyeLength;

        real cutOffDst;
            
        struct Parameters{

            real debyeLength;

            real cutOffDst;
        };

        DebyeHuckelDistanceDependentDielectric_(std::shared_ptr<ParticleGroup> pg,
                                                std::shared_ptr<Topology>     top,
                                                Parameters par):pg(pg),
                                                                pd(pg->getParticleData()),
                                                                sys(pg->getParticleData()->getSystem()),
                                                                debyeLength(par.debyeLength),
                                                                cutOffDst(par.cutOffDst){}
        
        real getCutOffDst(){
            return cutOffDst;
        }

        ~DebyeHuckelDistanceDependentDielectric_(){}

        struct forceTransverser{

            real4* force;
            
            real* charge;

            real debyeLength;

            Box box;

            real cutOffDst2;

            forceTransverser(real4* force,
                             real* charge,
                             real debyeLength,
                             Box box,
                             real cutOffDst2):force(force),
                                           charge(charge),
                                           debyeLength(debyeLength),
                                           box(box),
                                           cutOffDst2(cutOffDst2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                if(r2>0 and r2<=cutOffDst2){

                    const real chgProduct = charge[index_i]*charge[index_j];

                    return make_real4(CommonPotentials::DebyeHuckel::DebyeHuckelDistanceDependentDielectric::force<Topology::Units>(rij,r2,chgProduct,debyeLength,
                                                                                                                                    real(-8.5525),real(78.4)-real(-8.5525),
                                                                                                                                    real(7.7839),
                                                                                                                                    real(0.003627)),
                                                                                                                                    real(0.0));
                } 

                return make_real4(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     

            return forceTransverser(force,
                                    charge,
                                    debyeLength,
                                    box,
                                    cutOffDst*cutOffDst);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* charge;

            real debyeLength;

            Box box;

            real cutOffDst2;

            energyTransverser(real* energy,
                              real* charge,
                              real debyeLength,
                              Box box,
                              real cutOffDst2):energy(energy),
                                            charge(charge),
                                            debyeLength(debyeLength),
                                            box(box),
                                            cutOffDst2(cutOffDst2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=cutOffDst2){
                    
                    const real chgProduct = charge[index_i]*charge[index_j];

                    return CommonPotentials::DebyeHuckel::DebyeHuckelDistanceDependentDielectric::energy<Topology::Units>(rij,r2,chgProduct,debyeLength,
                                                                                                                          real(-8.5525),real(78.4)-real(-8.5525),
                                                                                                                          real(7.7839),
                                                                                                                          real(0.003627));
                } 

                return real(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy  = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
            
            return energyTransverser(energy,
                                     charge,
                                     debyeLength,
                                     box,
                                     cutOffDst*cutOffDst);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };
    
    template<typename Topology>
    using DebyeHuckelDistanceDependentDielectric = DebyeHuckelDistanceDependentDielectric_<Topology>;

}}}}

#endif
