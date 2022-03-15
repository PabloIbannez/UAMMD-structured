#ifndef __DH_UNBOUND__
#define __DH_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{

    template<typename Units_>
    struct DebyeHuckel: public ParameterUpdatable{
        
        std::shared_ptr<ParticleData> pd;

        Box box;
            
        real dielectricConstant; 
        real debyeLength;

        real cutOff;
            
        struct Parameters{

            real dielectricConstant; 
            real debyeLength;

            real cutOff;
        };

        DebyeHuckel(std::shared_ptr<ParticleData> pd,
                    Parameters par):pd(pd),
                                    dielectricConstant(par.dielectricConstant),
                                    debyeLength(par.debyeLength),
                                    cutOff(par.cutOff){}

        ~DebyeHuckel(){}

        struct forceTransverser{

            real4* force;
            
            real* charge;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOff2;

            forceTransverser(real4* force,
                             real* charge,
                             real dielectricConstant, 
                             real debyeLength,
                             Box box,
                             real cutOff2):force(force),
                                           charge(charge),
                                           dielectricConstant(dielectricConstant),
                                           debyeLength(debyeLength),
                                           box(box),
                                           cutOff2(cutOff2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);

                if(r2>0 and r2<=cutOff2){

                    const real chgProduct = charge[index_i]*charge[index_j];

                    return make_real4(CommonPotentials::DebyeHuckel::DebyeHuckel::force<Units_>(rij,r2,chgProduct,
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
                                    cutOff*cutOff);
        }
        
        struct virialTransverser{

            tensor3* virial;
                              
            real* charge;


            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOff2;

            virialTransverser(tensor3* virial,
                              real* charge,
                              real dielectricConstant, 
                              real debyeLength,
                              Box box,
                              real cutOff2):virial(virial),
                                            charge(charge),
                                            dielectricConstant(dielectricConstant),
                                            debyeLength(debyeLength),
                                            box(box),
                                            cutOff2(cutOff2){}

            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                    
                if(r2>0 and r2<=cutOff2){
                    
                    const real chgProduct = charge[index_i]*charge[index_j];
                    
                    return (CommonPotentials::DebyeHuckel::DebyeHuckel::virial<Units_>(rij,r2,chgProduct,
                                                                                       dielectricConstant,debyeLength));
                } 

                return tensor3(0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){virial[index_i]+=quantity;}

        };

        virialTransverser getVirialTransverser(){
            
            tensor3* virial = this->pd->getVirial(access::location::gpu, access::mode::readwrite).raw();     
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     

            return virialTransverser(virial,
                                     charge,
                                     dielectricConstant, 
                                     debyeLength,
                                     box,
                                     cutOff*cutOff);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* charge;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOff2;

            energyTransverser(real* energy,
                              real* charge,
                              real dielectricConstant, 
                              real debyeLength,
                              Box box,
                              real cutOff2):energy(energy),
                                            charge(charge),
                                            dielectricConstant(dielectricConstant),
                                            debyeLength(debyeLength),
                                            box(box),
                                            cutOff2(cutOff2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                if(r2>0 and r2<=cutOff2){
                    
                    const real chgProduct = charge[index_i]*charge[index_j];

                    return (CommonPotentials::DebyeHuckel::DebyeHuckel::energy<Units_>(rij,r2,chgProduct,
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
                                     cutOff*cutOff);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };

}}}}

#endif
