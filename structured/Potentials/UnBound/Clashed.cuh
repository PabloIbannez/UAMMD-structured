#ifndef __CLASHED_UNBOUND__
#define __CLASHED_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
    template<class Topology>
    struct Clashed_: public ParameterUpdatable{
        
        std::shared_ptr<System>        sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<Topology>     top;

        Box box;

        real lambda;
        real gamma;

        real cutOffDst;
            
        struct Parameters{
            
            real lambda;
            real gamma;

            real cutOffDst;
        
        };

        Clashed_(std::shared_ptr<ParticleGroup> pg,
                 std::shared_ptr<Topology>     top,
                 Parameters par):pg(pg),
                                 pd(pg->getParticleData()),
                                 sys(pg->getParticleData()->getSystem()),
                                 lambda(par.lambda),
                                 gamma(par.gamma),
                                 cutOffDst(par.cutOffDst){}

        real getCutOffDst(){
            return cutOffDst;
        }

        ~Clashed_(){}

        struct forceTransverser{

            real4* force;
            real* radius;

            real lambda;
            real gamma;
            
            Box  box;

            real cutOffDst2;
            
            forceTransverser(real4* force,
                             real* radius,
                             real lambda,
                             real gamma,
                             Box  box,
                             real cutOffDst2):force(force),
                                              radius(radius),
                                              lambda(lambda),
                                              gamma(gamma),
                                              box(box),
                                              cutOffDst2(cutOffDst2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real d = gamma*(radius[index_i]+radius[index_j]);
                
                real fmod = 0;

                if(r2>0 and r2<=cutOffDst2){
                    real d2 = d*d;

                    if(d2 > r2){
                        fmod = -real(4.0)*(d2-r2); //fmod*rij/r=-real(4.0)*r*(d2-r2)*rij/r
                    } else {
                        fmod = 0;
                    }
                }

                return make_real4(lambda*fmod*rij,real(0.0));
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            real*  radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();

            return forceTransverser(force,
                                    radius,
                                    lambda,
                                    gamma,
                                    box,
                                    cutOffDst*cutOffDst);
        }
        
        struct energyTransverser{

            real* energy;
            real* radius;

            real lambda;
            real gamma;
            
            Box  box;
            
            real cutOffDst2;
            
            energyTransverser(real* energy,
                              real* radius,
                              real lambda,
                              real gamma,
                              Box  box,
                              real cutOffDst2):energy(energy),
                                               radius(radius),
                                               lambda(lambda),
                                               gamma(gamma),
                                               box(box),
                                               cutOffDst2(cutOffDst2){}


            using resultType=real;

            inline __device__ resultType zero(){return real(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real d = gamma*(radius[index_i]+radius[index_j]);

                real e = real(0.0);

                if(r2>0 and r2<=cutOffDst2){
                    e = max(real(0.0),d*d-r2);
                    e = e*e;
                }

                return lambda*e/real(2.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            
            return energyTransverser(energy,
                                     radius,
                                     lambda,
                                     gamma,
                                     box,
                                     cutOffDst*cutOffDst);
        }
        
        void updateBox(Box newBox) override {
            box=newBox;
        }

    };

    template<class Topology>
    using Clashed = Clashed_<Topology>;

}}}}

#endif
