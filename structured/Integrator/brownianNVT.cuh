#ifndef __BROWNIAN__
#define __BROWNIAN__

namespace uammd{
namespace structured{

    namespace BrownianNVT{

        class Newton : public IntegratorBasicNVT{
            
            private:

                int N;

            public:

                struct Parameters: public IntegratorBasicNVT::Parameters{
                    real frictionConstant;
                };

            private:
                
                real frictionConstant;
            
            public:
                
                
                Parameters inputFileToParam(InputFile& in){
                    
                    Parameters param;
                    static_cast<IntegratorBasicNVT::Parameters&>(param) = IntegratorBasicNVT::inputFileToParam(in); 
                
                    in.getOption("frictionConstant",InputFile::Required)
                                  >>param.frictionConstant;

                    return param;
                }

                Newton(shared_ptr<ParticleGroup> pg,
                       uammd::InputFile& in,
                       cudaStream_t stream):Newton(pg,inputFileToParam(in),stream){}

                Newton(shared_ptr<ParticleGroup> pg,
                       Parameters param,
                       cudaStream_t stream):IntegratorBasicNVT(pg,param,"BrownianNVT::Newton",stream),
                                            frictionConstant(param.frictionConstant){

                            sys->log<System::MESSAGE>("[%s] frictionConstant: %f",this->name.c_str(),frictionConstant);

                            IntegratorBasic_ns::loadFrictionConstant(pg,frictionConstant);
                }

                real getFrictionConstant(){
                    return frictionConstant;
                }
                
                ~Newton(){}
                
                template<class UNITS>
                void applyUnits(){
                    IntegratorBasicNVT::applyUnits<UNITS>();

                    frictionConstant = frictionConstant/UNITS::TO_INTERNAL_TIME;
                    IntegratorBasic_ns::loadFrictionConstant(pg,frictionConstant);
                    
                    sys->log<System::MESSAGE>("[%s] FrictionConstant (after units): %f", this->name.c_str(), frictionConstant);
                    
                }
                
                void init(){
                    CudaCheckError();
                    this->update();
                    
                    sys->log<System::MESSAGE>("[%s] Performing initialization step", name.c_str());
                    
                    N = pg->getNumberParticles();

                }

                void forwardTime() {
                    
                    if(steps==0){
                        this->init();
                    }
                    
                    CudaCheckError();
                    this->update();

                    steps++;
                    sys->log<System::DEBUG1>("[%s] Performing integration step %d", name.c_str(), steps);
                    
                    this->updateForce();
                    this->integrationStep();
                }
                        
                void integrationStep(){
                    
                    auto mass     = pd->getMass(access::location::gpu, access::mode::read);
                    auto friConst = pd->getFrictionConstant(access::location::gpu, access::mode::read);     
                    
                    auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                    
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);
                    
                    real* mass_ptr = mass.raw();
                    real* friConst_ptr = friConst.raw();
                    
                    real4* pos_ptr   = pos.raw();
                    real4* force_ptr = force.raw();
                
                    uint step_temp = steps;
                    uint seed_temp = seed;
                    real dt_temp = dt;
                    real kBT_temp = kB*T;
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){real fc = friConst_ptr[index]*mass_ptr[index]; // gamma_i = gamma_i*m_i

                                                                real sigma = sqrt(real(2.0)*kBT_temp*fc/dt_temp);

                                                                Saru rng_s(index, step_temp, seed_temp);
                                                                const real3 noise = make_real3(rng_s.gf(0, real(1.0)), rng_s.gf(0,real(1.0)).x);

                                                                real3 f = (dt_temp/fc)*(make_real3(force_ptr[index])+sigma*noise);

                                                                pos_ptr[index].x+=f.x; //pos also stores type at .w
                                                                pos_ptr[index].y+=f.y;
                                                                pos_ptr[index].z+=f.z;
                                                                
                                                                force_ptr[index] = make_real4(0);});
                }
        };
        
        class EulerMaruyama : public IntegratorBasicNVT{
            
            private:

                int N;

            public:

                struct Parameters: public IntegratorBasicNVT::Parameters{
                    real viscosity;
                };

            private:
                
                real viscosity;
            
            public:
                
                
                Parameters inputFileToParam(InputFile& in){
                    
                    Parameters param;
                    static_cast<IntegratorBasicNVT::Parameters&>(param) = IntegratorBasicNVT::inputFileToParam(in); 
                
                    in.getOption("viscosity",InputFile::Required)
                                  >>param.viscosity;

                    return param;
                }

                EulerMaruyama(shared_ptr<ParticleGroup> pg,
                              uammd::InputFile& in,
                              cudaStream_t stream):EulerMaruyama(pg,inputFileToParam(in),stream){}

                EulerMaruyama(shared_ptr<ParticleGroup> pg,
                              Parameters param,
                              cudaStream_t stream):IntegratorBasicNVT(pg,param,"BrownianNVT::EulerMaruyama",stream),
                                                   viscosity(param.viscosity){
  
                              sys->log<System::MESSAGE>("[%s] viscosity: %f",this->name.c_str(),viscosity);
  
                              loadMobility();
                           }
                
                ~EulerMaruyama(){}

                void loadMobility(){
                    
                    auto radius = pd->getRadius(access::location::cpu, access::mode::read);
                                  
                    auto mt = pd->getTranslationalSelfDiffusion(access::location::cpu, access::mode::write);
                    auto mr = pd->getRotationalSelfDiffusion(access::location::cpu, access::mode::write);
            
                    auto id = pd->getId(access::location::cpu, 
                                       access::mode::read);
            
                    auto groupIndex  = pg->getIndexIterator(access::location::cpu);
                    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

                    fori(0,pg->getNumberParticles()){
                        int index = sortedIndex[id[groupIndex[i]]];
                        mt[index] = real(1.0)/(real(6.0)*M_PI*viscosity*radius[index]);
                        mr[index] = real(1.0)/(real(8.0)*M_PI*viscosity*radius[index]*radius[index]*radius[index]);
                    }

                }
                
                template<class UNITS>
                void applyUnits(){
                    IntegratorBasicNVT::applyUnits<UNITS>();
                    
                    viscosity = viscosity/UNITS::TO_INTERNAL_TIME; //Not checked!! 
                    
                    loadMobility();
                    sys->log<System::MESSAGE>("[%s] Viscosity (after units): %f", this->name.c_str(), viscosity);
                    
                }
                
                void init(){
                    CudaCheckError();
                    this->update();
                    
                    sys->log<System::MESSAGE>("[%s] Performing initialization step", name.c_str());
                    
                    N = pg->getNumberParticles();
                }

                void forwardTime() {
                    
                    if(steps==0){
                        this->init();
                    }
                    
                    CudaCheckError();
                    this->update();

                    steps++;
                    sys->log<System::DEBUG1>("[%s] Performing integration step %d", name.c_str(), steps);
                    
                    this->updateForce();
                    this->integrationStep();
                }
                
                void integrationStep(){
                    
                    auto mt = pd->getTranslationalSelfDiffusion(access::location::gpu, access::mode::read);
                    auto mr = pd->getRotationalSelfDiffusion(access::location::gpu, access::mode::read);
                    
                    auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                    
                    auto dir    = pd->getDir(access::location::gpu, access::mode::readwrite);     
                    auto torque = pd->getTorque(access::location::gpu, access::mode::readwrite);     
                    
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);
                    
                    real* mt_ptr = mt.raw();
                    real* mr_ptr = mr.raw();
                    
                    real4* pos_ptr   = pos.raw();
                    real4* force_ptr = force.raw();
                    
                    real4* dir_ptr    = dir.raw();
                    real4* torque_ptr = torque.raw();

                    uint step_temp = steps;
                    uint seed_temp = seed;
                    real dt_temp = dt;
                    real kBT_temp = kB*T;

                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){
                                                                Saru rng_s(index, step_temp, seed_temp);

                                                                real2 rndAux = rng_s.gf(0, real(1.0));
                                                                real3 noiseTrans = make_real3(rng_s.gf(0, real(1.0)),rndAux.x);
                                                                real3 noiseRot   = make_real3(rng_s.gf(0, real(1.0)),rndAux.y);
                                                                
                                                                //Translation
                                                                real sigma; 
                                                                sigma = sqrt(real(2.0)*kBT_temp*mt_ptr[index]/dt_temp);

                                                                real3 f = make_real3(force_ptr[index]);
                                                                real3 v = mt_ptr[index]*f+sigma*noiseTrans;
                                                                
                                                                pos_ptr[index].x+=dt_temp*v.x;
                                                                pos_ptr[index].y+=dt_temp*v.y;
                                                                pos_ptr[index].z+=dt_temp*v.z;

                                                                //Rotation
                                                                sigma = sqrt(real(2.0)*kBT_temp*mr_ptr[index]/dt_temp);

                                                                real3 torque = make_real3(torque_ptr[index]);
                                                                real3 domega=mr_ptr[index]*torque+sigma*noiseRot;

                                                                dir[index] = quaternions::rotate(dir[index],dt_temp*domega);
                                                                dir[index] = normalize(dir[index]);
                                                                
                                                                force_ptr[index]  = make_real4(0);
                                                                torque_ptr[index] = make_real4(0);});
                }
        };
    }
    
}}

#endif

