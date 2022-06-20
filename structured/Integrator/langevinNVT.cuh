#ifndef __LANGEVIN_NVT__
#define __LANGEVIN_NVT__

namespace uammd{
namespace structured{
                
    namespace LangevinNVT{
    
    class BBK: public IntegratorBasicNVT{
        
        private:
            
            //System variables

            int N;
            real totalMass;

            //pos ref and vel ref stores the position and velocity of the previous integration step.
            //These states are needed, for example, to compute the temperature, for which you need the velocity at t+dt,
            //the velocity at t+dt can be computed used the velocity at t+3/2dt (stored at pd->vel) an the velocity at t+1/2dt (stored at velRef).
            //For the stopTransRot algorithm you have to use the positon at t+0, which stored at posRef
            thrust::device_vector<real4> posRef;
            thrust::device_vector<real3> velRef;
            
            thrust::device_vector<real4> posTemp;

        public:

            struct Parameters: public IntegratorBasicNVT::Parameters{
                real frictionConstant;
            };

        private:
            
            real frictionConstant;
            
            void copyToRef(){
                
                auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);     
                posRef.resize(pd->getNumParticles());
                thrust::copy(thrust::cuda::par.on(stream),pos.begin(), pos.end(), posRef.begin());
                
                auto vel = pd->getVel(access::location::gpu, access::mode::readwrite);     
                velRef.resize(pd->getNumParticles());
                thrust::copy(thrust::cuda::par.on(stream),vel.begin(), vel.end(), velRef.begin());

            }
            
            void copyFromRef(){
                
                auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);     
                thrust::copy(thrust::cuda::par.on(stream),posRef.begin(), posRef.end(), pos.begin());
                
                auto vel = pd->getVel(access::location::gpu, access::mode::readwrite);     
                thrust::copy(thrust::cuda::par.on(stream),velRef.begin(), velRef.end(), vel.begin());

            }
            
            Parameters inputFileToParam(InputFile& in){
                
                Parameters param;
                static_cast<IntegratorBasicNVT::Parameters&>(param) = IntegratorBasicNVT::inputFileToParam(in); 
            
                in.getOption("frictionConstant",InputFile::Required)
                              >>param.frictionConstant;

                return param;
            }

        public:
            
            BBK(shared_ptr<ParticleGroup> pg,
                uammd::InputFile& in,
                cudaStream_t stream):BBK(pg,inputFileToParam(in),stream){}
            
            BBK(shared_ptr<ParticleGroup> pg,
                Parameters param,
                cudaStream_t stream):IntegratorBasicNVT(pg,param, "LangevinNVT::BBK",stream),
                                     frictionConstant(param.frictionConstant){

                sys->log<System::MESSAGE>("[%s] frictionConstant: %f",this->name.c_str(),frictionConstant);

                IntegratorBasic_ns::loadFrictionConstant(pg,frictionConstant);

            }

            real getFrictionConstant(){
                return frictionConstant;
            }

            template<class UNITS>
            void applyUnits(){
                IntegratorBasicNVT::applyUnits<UNITS>();
                
                frictionConstant = frictionConstant/UNITS::TO_INTERNAL_TIME;
                IntegratorBasic_ns::loadFrictionConstant(pg,frictionConstant);
                
                sys->log<System::MESSAGE>("[%s] FrictionConstant (after units): %f", this->name.c_str(), frictionConstant);
                
            }

            void stopTransRot(){
                
                //Remove translational and rotational motion about COM at (t + 1/2dt)

                sys->log<System::DEBUG1>("[%s] Stopping translation and rotation ", name.c_str());
                //sys->log<System::MESSAGE>("[%s] Stopping translation and rotation ", name.c_str());

                posTemp.resize(pg->getNumberParticles());
                
                {
                    auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
                    thrust::copy(thrust::cuda::par.on(stream),pos.begin(), pos.end(), posTemp.begin());
                    //CudaSafeCall(cudaStreamSynchronize(stream));
                }

                {
                    //Compute 0.5(posTemp+posRef) and store it in particle data position
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);
                    
                    auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);

                    real4* pos_ptr     = pos.raw();
                    real4* posTemp_ptr = thrust::raw_pointer_cast(posTemp.data());
                    real4* posRef_ptr  = thrust::raw_pointer_cast(posRef.data());
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){pos_ptr[index]=real(0.5)*(posTemp_ptr[index]+posRef_ptr[index]);});
                    //CudaSafeCall(cudaStreamSynchronize(stream));
                    //pos      stores pos(t+0.5dt)
                    //posTemp  stores pos(t+dt)
                    //posRef   stores pos(t+0)
                }

                //Compute the center of mass postion and velocity
                real3 comp = Measures::centerOfMassPos(pg,totalMass,stream);
                real3 comv = Measures::centerOfMassVel(pg,totalMass,stream);
                
                //Compute the center of mass angular velocity
                real3 angv = Measures::angularVelocity(pg,comp,comv,stream);
                
                {
                    
                    //Substract global movements
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);
                    
                    auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
                    auto vel = pd->getVel(access::location::gpu, access::mode::readwrite);
                    
                    real4* pos_ptr=pos.raw();
                    real3* vel_ptr=vel.raw();
                    
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){vel_ptr[index]=vel_ptr[index]-comv;
                                                                real3 d_comp_part = make_real3(pos_ptr[index])-comp;
                                                                vel_ptr[index].x=vel_ptr[index].x - (angv.y*d_comp_part.z-angv.z*d_comp_part.y);
                                                                vel_ptr[index].y=vel_ptr[index].y - (angv.z*d_comp_part.x-angv.x*d_comp_part.z);
                                                                vel_ptr[index].z=vel_ptr[index].z - (angv.x*d_comp_part.y-angv.y*d_comp_part.x);
                            });
                    
                    //CudaSafeCall(cudaStreamSynchronize(stream));
                }
                
                {
                    //Copy back posTemp to pos, now pos stores pos(t+dt)
                    auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
                    thrust::copy(thrust::cuda::par.on(stream),posTemp.begin(), posTemp.end(), pos.begin());
                    //CudaSafeCall(cudaStreamSynchronize(stream));
                }
            }

            void resetVelocities(){
                IntegratorBasic_ns::generateVelocity(this->pg,
                                                     this->kB*this->T,
                                                     this->name,
                                                     this->stream);
            }
            
            void init(){
                CudaCheckError();
                this->update();
                
                sys->log<System::MESSAGE>("[%s] Performing initialization step", name.c_str());
                
                N = pg->getNumberParticles();
                
                totalMass = Measures::totalMass(pg,stream);
                
                //refPos(0) refVel(0)
                this->copyToRef();
                
                //Inital velocities for some T
                this->resetVelocities();
                
                if(stopTransRotSteps > 0){
                    //Stop global movements
                    this->stopTransRot();
                }
                
                this->copyToRef();
                
                if(constrained){
                    this->constraint->applyConstraint(stream);
                    this->copyToRef();
                }
                
                //Reset force to 0
                this->resetForce();
                this->updateForce();
                //CudaSafeCall(cudaStreamSynchronize(stream));
                
                if(constrained){
                 
                    {
                        auto mass = pd->getMass(access::location::gpu, access::mode::read);
                        
                        auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                        auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     

                        auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                        
                        auto groupIterator = pg->getIndexIterator(access::location::gpu);
                        
                        real* mass_ptr    = mass.raw();
                        
                        real4* pos_ptr    = pos.raw();
                        real4* posRef_ptr = thrust::raw_pointer_cast(posRef.data());
                        real3* vel_ptr    = vel.raw();
                        real3* velRef_ptr = thrust::raw_pointer_cast(velRef.data());

                        real4* force_ptr  = force.raw();
                        
                        real dt_temp = dt;
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){const real invMass = real(1.0)/mass_ptr[index];
                                                                    vel_ptr[index] = velRef_ptr[index]+real(0.5)*dt_temp*make_real3(force_ptr[index])*invMass;

                                                                    pos_ptr[index].x = posRef_ptr[index].x+dt_temp*vel_ptr[index].x;
                                                                    pos_ptr[index].y = posRef_ptr[index].y+dt_temp*vel_ptr[index].y;
                                                                    pos_ptr[index].z = posRef_ptr[index].z+dt_temp*vel_ptr[index].z;
                                                                    
                                                                    force_ptr[index]=make_real4(0.0);});
                    }
                    
                    this->constraint->applyConstraint(stream);

                    {
                        auto mass = pd->getMass(access::location::gpu, access::mode::read);
                        
                        auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                        auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     

                        auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                        
                        auto groupIterator = pg->getIndexIterator(access::location::gpu);
                        
                        real* mass_ptr    = mass.raw();
                        
                        real4* pos_ptr    = pos.raw();
                        real4* posRef_ptr = thrust::raw_pointer_cast(posRef.data());
                        real3* vel_ptr    = vel.raw();
                        real3* velRef_ptr = thrust::raw_pointer_cast(velRef.data());

                        real4* force_ptr  = force.raw();
                    
                        real dt_temp = dt;
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){
                                                                    vel_ptr[index]=(make_real3(pos_ptr[index])-make_real3(posRef_ptr[index]))/dt_temp;
                                                                    });
                        
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){const real invMass = real(1.0)/mass_ptr[index];
                                                                    velRef_ptr[index] = velRef_ptr[index]-real(0.5)*dt_temp*make_real3(force_ptr[index])*invMass;

                                                                    pos_ptr[index].x = posRef_ptr[index].x-dt_temp*vel_ptr[index].x;
                                                                    pos_ptr[index].y = posRef_ptr[index].y-dt_temp*vel_ptr[index].y;
                                                                    pos_ptr[index].z = posRef_ptr[index].z-dt_temp*vel_ptr[index].z;
                                                                    
                                                                    force_ptr[index]=make_real4(0.0);});
                    }
                    
                    this->constraint->applyConstraint(stream);
                       
                    {
                        auto mass = pd->getMass(access::location::gpu, access::mode::read);
                        
                        auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                        auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     

                        auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                        
                        auto groupIterator = pg->getIndexIterator(access::location::gpu);
                        
                        real* mass_ptr    = mass.raw();
                        
                        real4* pos_ptr    = pos.raw();
                        real4* posRef_ptr = thrust::raw_pointer_cast(posRef.data());
                        real3* vel_ptr    = vel.raw();
                        real3* velRef_ptr = thrust::raw_pointer_cast(velRef.data());

                        real4* force_ptr  = force.raw();
                    
                        real dt_temp = dt;
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){
                                                                    velRef_ptr[index]=(make_real3(posRef_ptr[index])-make_real3(pos_ptr[index]))/dt_temp;
                                                                    });
                    }

                    this->copyFromRef();

                } else {
                
                    auto mass = pd->getMass(access::location::gpu, access::mode::read);
                    
                    auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                    auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                    
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);
                    
                    real* mass_ptr    = mass.raw();
                    
                    real4* pos_ptr   = pos.raw();
                    real3* vel_ptr   = vel.raw();
                    real4* force_ptr = force.raw();
            
                    real dt_temp = dt;
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){const real invMass = real(1.0)/mass_ptr[index];
                                                                vel_ptr[index]  =vel_ptr[index]-real(0.5)*dt_temp*make_real3(force_ptr[index])*invMass;
                                                                force_ptr[index]=make_real4(0.0);});
                    
                    //CudaSafeCall(cudaStreamSynchronize(stream));
                }
                
                //refPos(0) refVel(0-dt/2)
                this->copyToRef();

                this->updateForce();

                CudaCheckError();
                {
                    auto mass = pd->getMass(access::location::gpu, access::mode::read);
                    
                    auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                    auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                    
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);
                    
                    real* mass_ptr    = mass.raw();
                    
                    real4* pos_ptr   = pos.raw();
                    real3* vel_ptr   = vel.raw();
                    real4* force_ptr = force.raw();
            
                    real dt_temp = dt;
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){const real invMass = real(1.0)/mass_ptr[index];

                                                                vel_ptr[index]=vel_ptr[index]+dt_temp*make_real3(force_ptr[index])*invMass;
                                                                
                                                                pos_ptr[index].x=pos_ptr[index].x+dt_temp*vel_ptr[index].x;
                                                                pos_ptr[index].y=pos_ptr[index].y+dt_temp*vel_ptr[index].y;
                                                                pos_ptr[index].z=pos_ptr[index].z+dt_temp*vel_ptr[index].z;

                                                                force_ptr[index]=make_real4(0.0);});
                    //CudaSafeCall(cudaStreamSynchronize(stream));
                }

                if(constrained){

                    this->constraint->applyConstraint(stream);

                    auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                    auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     
                    
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);
                    
                    real4* pos_ptr    = pos.raw();
                    real4* posRef_ptr = thrust::raw_pointer_cast(posRef.data());
                    real3* vel_ptr    = vel.raw();
            
                    real dt_temp = dt;
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){
                                                                vel_ptr[index]=(make_real3(pos_ptr[index])-make_real3(posRef_ptr[index]))/dt_temp;
                                                                });

                    //TODO virial constraint
                }
            }
            
            void forwardTime() {

                if(steps==0){
                    this->init();
                }
                
                CudaCheckError();
                this->update();

                steps++;
                sys->log<System::DEBUG1>("[%s] Performing integration step %d", name.c_str(), steps);
                
                
                if(stopTransRotSteps > 0 or constrained){
                    if(!constrained) {
                        if(steps%stopTransRotSteps == 0){
                            this->copyToRef();
                        }
                    } else {
                        this->copyToRef();
                    }
                }
                
                this->updateForce();
                this->integrationStep();

                if(stopTransRotSteps > 0){
                    if(steps%stopTransRotSteps == 0){
                        //Stop global movements
                        this->stopTransRot();
                    }
                }
            }
                    
            void integrationStep(){

                int numberParticles = pg->getNumberParticles();
                
                auto mass     = pd->getMass(access::location::gpu, access::mode::read);
                auto friConst = pd->getFrictionConstant(access::location::gpu, access::mode::read);     
                
                auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     
                auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                
                auto groupIterator = pg->getIndexIterator(access::location::gpu);
                
                real* mass_ptr = mass.raw();
                real* friConst_ptr = friConst.raw();
                
                real4* pos_ptr   = pos.raw();
                real3* vel_ptr   = vel.raw();
                real4* force_ptr = force.raw();
            
                uint step_temp = steps;
                uint seed_temp = seed;
                real dt_temp  = dt;
                real kBT_temp = kB*T;
                thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                        [=] __host__ __device__ (int index){const real fricConstant = friConst_ptr[index];
                                                            
                                                            real factor  = real(2.0)*fricConstant*kBT_temp/dt_temp;
                                                            real factor2 = real(1.0)+real(0.5)*dt_temp*fricConstant;
                                                            
                                                            const real m = mass_ptr[index];
                                                            
                                                            
                                                            Saru rng_s(index, step_temp, seed_temp);
                                                            const real3 noise = make_real3(rng_s.gf(0, real(1.0)), rng_s.gf(0,real(1.0)).x);
                
                                                            real scaleV  = real(2.0)/factor2 - real(1.0);
                                                            real scaleF  = dt_temp/factor2;
                                                            
                                                            real sigma   = sqrt(m*factor);
                                                            real3 f = make_real3(force_ptr[index])+sigma*noise;
                                                            
                                                            const real invMass = real(1.0)/m;
                                                            vel[index]=scaleV*vel[index]+scaleF*f*invMass;

                                                            pos[index].x=pos[index].x+dt_temp*vel[index].x; //pos also stores type at .w
                                                            pos[index].y=pos[index].y+dt_temp*vel[index].y;
                                                            pos[index].z=pos[index].z+dt_temp*vel[index].z;
                                                            
                                                            force[index] = make_real4(0);});
                
                if(constrained){

                    this->constraint->applyConstraint(stream);

                    real4* posRef_ptr = thrust::raw_pointer_cast(posRef.data());
                    
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){
                                                                vel_ptr[index]=(make_real3(pos_ptr[index])-make_real3(posRef_ptr[index]))/dt_temp;
                                                                });

                    //TODO virial constraint
                }

                //CudaSafeCall(cudaStreamSynchronize(stream));
                CudaCheckError();

            }
    };
    
    class GJF: public IntegratorBasicNVT{
        
        private:
            
            //System variables

            int N;

        public:

            struct Parameters: public IntegratorBasicNVT::Parameters{
                real frictionConstant;
            };

        private:
            
            real frictionConstant;
            
            Parameters inputFileToParam(InputFile& in){
                
                Parameters param;
                static_cast<IntegratorBasicNVT::Parameters&>(param) = IntegratorBasicNVT::inputFileToParam(in); 
            
                in.getOption("frictionConstant",InputFile::Required)
                              >>param.frictionConstant;

                return param;
            }

        public:
            
            GJF(shared_ptr<ParticleGroup> pg,
                uammd::InputFile& in,
                cudaStream_t stream):GJF(pg,inputFileToParam(in),stream){}
            
            GJF(shared_ptr<ParticleGroup> pg,
                Parameters param,
                cudaStream_t stream):IntegratorBasicNVT(pg, param, "LangevinNVT::GJF",stream),
                                     frictionConstant(param.frictionConstant){

                sys->log<System::MESSAGE>("[%s] frictionConstant: %f",this->name.c_str(),frictionConstant);

                IntegratorBasic_ns::loadFrictionConstant(pg,frictionConstant);

            }
            
            template<class UNITS>
            void applyUnits(){
                IntegratorBasicNVT::applyUnits<UNITS>();
                
                frictionConstant = frictionConstant/UNITS::TO_INTERNAL_TIME;
                IntegratorBasic_ns::loadFrictionConstant(pg,frictionConstant); 
                
                sys->log<System::MESSAGE>("[%s] FrictionConstant (after units): %f", this->name.c_str(), frictionConstant);
                
            }
            
            void resetVelocities(){
                IntegratorBasic_ns::generateVelocity(this->pg,
                                                     this->kB*this->T,
                                                     this->name,
                                                     this->stream);
            }
            
            void init(){
                CudaCheckError();
                this->update();
                
                sys->log<System::MESSAGE>("[%s] Performing initialization step", name.c_str());
                
                N = pg->getNumberParticles();
                
                //Inital velocities for some T
                this->resetVelocities();
                
                this->resetForce();
                this->updateForce();
                
            }
            
            void forwardTime() {

                if(steps==0){
                    this->init();
                }
                
                CudaCheckError();
                this->update();

                steps++;
                sys->log<System::DEBUG1>("[%s] Performing integration step %d", name.c_str(), steps);
                
                this->integrationStep();
            }
                    
            void integrationStep(){

                int numberParticles = pg->getNumberParticles();
                
                uint step_temp  = steps;
                uint seed_temp  = seed;
                                                            
                real comboGJ = 0.5*dt*frictionConstant;
	            real sigmaGaussPrefactor = sqrt(2*frictionConstant*kB*T*dt);
                real bGJ = 1.0/(1.0+comboGJ);
                real aGJ = (1.0-comboGJ)/(1.0+comboGJ);
	            real cr2 = bGJ*dt;
	            real cr3 = bGJ*dt*dt*0.5; 
	            real cr4 = bGJ*dt*0.5; 
	            real vr1 = aGJ; 
	            real vr2 = dt*0.5*aGJ;  
	            real vr3 = dt*0.5;  
	            real vr4 = bGJ;
                
                {
                    auto mass  = pd->getMass(access::location::gpu, access::mode::read);
                    
                    auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite);     
                    auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                    
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);

                    real*  mass_ptr = mass.raw();
                    
                    real4* pos_ptr   = pos.raw();
                    real3* vel_ptr   = vel.raw();
                    real4* force_ptr = force.raw();
        
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){                      
                                                                const real m = mass_ptr[index];

                                                                real sigma   = sqrt(m)*sigmaGaussPrefactor;
                                                                
                                                                Saru rng_s(index, step_temp, seed_temp);
                                                                real3 noise = make_real3(rng_s.gf(real(0.0), real(1.0)), rng_s.gf(real(0.0),real(1.0)).x);

                                                                noise = sigma*noise;

                                                                real3 f = make_real3(force_ptr[index]);

                                                                pos_ptr[index].x=pos_ptr[index].x+cr2*vel_ptr[index].x+(cr3/m)*f.x+(cr4/m)*noise.x; //pos also stores type at .w
                                                                pos_ptr[index].y=pos_ptr[index].y+cr2*vel_ptr[index].y+(cr3/m)*f.y+(cr4/m)*noise.y;
                                                                pos_ptr[index].z=pos_ptr[index].z+cr2*vel_ptr[index].z+(cr3/m)*f.z+(cr4/m)*noise.z;
                                                                
                                                                vel_ptr[index].x=vr1*vel_ptr[index].x+(vr2/m)*f.x+(vr4/m)*noise.x;
                                                                vel_ptr[index].y=vr1*vel_ptr[index].y+(vr2/m)*f.y+(vr4/m)*noise.y;
                                                                vel_ptr[index].z=vr1*vel_ptr[index].z+(vr2/m)*f.z+(vr4/m)*noise.z;
                                                                
                                                                force_ptr[index] = make_real4(0);});
                }

                this->updateForce();
                
                {
                    auto mass  = pd->getMass(access::location::gpu, access::mode::read);
                    
                    auto vel   = pd->getVel(access::location::gpu, access::mode::readwrite);     
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                    
                    auto groupIterator = pg->getIndexIterator(access::location::gpu);

                    real*  mass_ptr = mass.raw();
                    
                    real3* vel_ptr   = vel.raw();
                    real4* force_ptr = force.raw();
                
                    thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                            [=] __host__ __device__ (int index){                      
                                                                const real m = mass_ptr[index];

                                                                real3 f = make_real3(force_ptr[index]);

                                                                vel_ptr[index].x=vel_ptr[index].x+(vr3/m)*f.x;
                                                                vel_ptr[index].y=vel_ptr[index].y+(vr3/m)*f.y;
                                                                vel_ptr[index].z=vel_ptr[index].z+(vr3/m)*f.z;
                                                                
                                                                });
                }

            }
            
    };

    }

}}

#endif

