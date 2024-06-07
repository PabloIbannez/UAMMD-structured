#ifndef __BBK_NVT__
#define __BBK_NVT__

namespace uammd{
namespace structured{
namespace NVT{

    namespace Langevin{

        class BBK: public IntegratorBasicNVT{

            private:

                ////////////////////////////////////

                ullint stopTransRotSteps;

                ////////////////////////////////////

                bool initialized     = false;
                bool resetVelocities = true;

                real totalMass;

                //pos ref and vel ref stores the position and velocity of the previous integration step.
                //These states are needed, for example, to compute the temperature, for which you need the velocity at t+dt,
                //the velocity at t+dt can be computed used the velocity at t+3/2dt (stored at pd->vel) an the velocity at t+1/2dt (stored at velRef).
                //For the stopTransRot algorithm you have to use the positon at t+0, which stored at posRef
                thrust::device_vector<real4> posRef;
                thrust::device_vector<real3> velRef;

                thrust::device_vector<real4> posTemp;

                ////////////////////////////////////


                void copyToRef(){

                    auto pos = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                    posRef.resize(this->pd->getNumParticles());
                    thrust::copy(thrust::cuda::par.on(stream),pos.begin(), pos.end(), posRef.begin());

                    auto vel = this->pd->getVel(access::location::gpu, access::mode::readwrite);
                    velRef.resize(this->pd->getNumParticles());
                    thrust::copy(thrust::cuda::par.on(stream),vel.begin(), vel.end(), velRef.begin());

                    CudaSafeCall(cudaStreamSynchronize(stream));

                }

                void copyFromRef(){

                    auto pos = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                    thrust::copy(thrust::cuda::par.on(stream),posRef.begin(), posRef.end(), pos.begin());

                    auto vel = this->pd->getVel(access::location::gpu, access::mode::readwrite);
                    thrust::copy(thrust::cuda::par.on(stream),velRef.begin(), velRef.end(), vel.begin());

                    CudaSafeCall(cudaStreamSynchronize(stream));

                }

            public:

                BBK(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data,
                    std::string name):IntegratorBasicNVT(gd,pg,data,name){

                    System::log<System::MESSAGE>("[BBK] Created BBK integrator \"%s\"",name.c_str());

                    if(!data.isParameterAdded("frictionConstant")){
                        if(this->pd->isFrictionConstantAllocated()){
                            System::log<System::WARNING>("[BBK] (%s) No friction constant specified, using the one in the particle data",name.c_str());
                        } else {
                            System::log<System::CRITICAL>("[BBK] (%s) No friction constant specified and none in the particle data",name.c_str());
                        }
                    } else {
                        real frictionConstant = data.getParameter<real>("frictionConstant");
                        System::log<System::MESSAGE>("[BBK] (%s) frictionConstant : %f",name.c_str(),frictionConstant);
                        IntegratorUtils::loadFrictionConstant(pg,frictionConstant);
                    }

                    stopTransRotSteps = data.getParameter<ullint>("stopTransRotSteps",0);
                    resetVelocities   = data.getParameter<bool>("resetVelocities",true);

                    System::log<System::MESSAGE>("[BBK] (%s) stopTransRotSteps : %llu",name.c_str(),stopTransRotSteps);
                    System::log<System::MESSAGE>("[BBK] (%s) resetVelocities : %i",name.c_str(),resetVelocities);
                }

                void stopTransRot(){

                    //Remove translational and rotational motion about COM at (t + 1/2dt)

                    System::log<System::DEBUG1>("[BBK] Stopping translation and rotation ");

                    int N = this->pg->getNumberParticles();

                    posTemp.resize(N);

                    {
                        auto pos = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                        thrust::copy(thrust::cuda::par.on(stream),pos.begin(), pos.end(), posTemp.begin());
                        CudaSafeCall(cudaStreamSynchronize(stream));
                    }

                    {
                        //Compute 0.5(posTemp+posRef) and store it in particle data position
                        auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                        auto pos = this->pd->getPos(access::location::gpu, access::mode::readwrite);

                        real4* pos_ptr     = pos.raw();
                        real4* posTemp_ptr = thrust::raw_pointer_cast(posTemp.data());
                        real4* posRef_ptr  = thrust::raw_pointer_cast(posRef.data());
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){pos_ptr[index]=real(0.5)*(posTemp_ptr[index]+posRef_ptr[index]);});
                        //pos      stores pos(t+0.5dt)
                        //posTemp  stores pos(t+dt)
                        //posRef   stores pos(t+0)
                        CudaSafeCall(cudaStreamSynchronize(stream));
                    }

                    //Compute the center of mass postion and velocity
                    real3 comp = Measures::centerOfMassPos(this->pg,totalMass,stream);
                    real3 comv = Measures::centerOfMassVel(this->pg,totalMass,stream);

                    //Compute the center of mass angular velocity
                    Box box = this->gd->getEnsemble()->getBox();
                    real3 angv = Measures::angularVelocity(this->pg,comp,comv,box,stream);

                    {

                        //Substract global movements
                        auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                        auto pos = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                        auto vel = this->pd->getVel(access::location::gpu, access::mode::readwrite);

                        real4* pos_ptr=pos.raw();
                        real3* vel_ptr=vel.raw();

                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){vel_ptr[index]=vel_ptr[index]-comv;
                                                                    real3 d_comp_part = make_real3(pos_ptr[index])-comp;
                                                                    vel_ptr[index].x=vel_ptr[index].x - (angv.y*d_comp_part.z-angv.z*d_comp_part.y);
                                                                    vel_ptr[index].y=vel_ptr[index].y - (angv.z*d_comp_part.x-angv.x*d_comp_part.z);
                                                                    vel_ptr[index].z=vel_ptr[index].z - (angv.x*d_comp_part.y-angv.y*d_comp_part.x);
                                });

                        CudaSafeCall(cudaStreamSynchronize(stream));
                    }

                    {
                        //Copy back posTemp to pos, now pos stores pos(t+dt)
                        auto pos = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                        thrust::copy(thrust::cuda::par.on(stream),posTemp.begin(), posTemp.end(), pos.begin());
                        CudaSafeCall(cudaStreamSynchronize(stream));
                    }

                    CudaCheckError();
                }

                void init(){

                    System::log<System::MESSAGE>("[BBK] Performing initialization step");

                    int N = this->pg->getNumberParticles();

                    totalMass = Measures::totalMass(this->pg);

                    //refPos(0) refVel(0)
                    this->copyToRef();

                    //Inital velocities for some T
                    if(resetVelocities){
                        IntegratorUtils::generateVelocity(this->pg,
                                                          this->kBT,
                                                          this->gd->getSystem()->getSeed(),
                                                          this->stream);
                    }

                    if(stopTransRotSteps > 0){
                        //Stop global movements
                        this->stopTransRot();
                    }

                    this->copyToRef();

                    //Reset force to 0
                    this->resetForce();
                    this->updateForce();
                    //CudaSafeCall(cudaStreamSynchronize(stream));
                    {
                        auto mass = this->pd->getMass(access::location::gpu, access::mode::read);

                        auto pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                        auto vel   = this->pd->getVel(access::location::gpu, access::mode::readwrite);
                        auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                        auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                        real* mass_ptr    = mass.raw();

                        real4* pos_ptr   = pos.raw();
                        real3* vel_ptr   = vel.raw();
                        real4* force_ptr = force.raw();

                        real dt_temp = this->dt;
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){const real invMass = real(1.0)/mass_ptr[index];
                                                                    vel_ptr[index]  =vel_ptr[index]-real(0.5)*dt_temp*make_real3(force_ptr[index])*invMass;
                                                                    force_ptr[index]=make_real4(0.0);});
                    }

                    CudaSafeCall(cudaStreamSynchronize(stream));

                    ////refPos(0) refVel(0-dt/2)
                    this->copyToRef();

                    this->updateForce();

                    CudaCheckError();
                    {
                        auto mass = this->pd->getMass(access::location::gpu, access::mode::read);

                        auto pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                        auto vel   = this->pd->getVel(access::location::gpu, access::mode::readwrite);
                        auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                        auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                        real* mass_ptr    = mass.raw();

                        real4* pos_ptr   = pos.raw();
                        real3* vel_ptr   = vel.raw();
                        real4* force_ptr = force.raw();

                        real dt_temp = this->dt;
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){const real invMass = real(1.0)/mass_ptr[index];

                                                                    vel_ptr[index]=vel_ptr[index]+dt_temp*make_real3(force_ptr[index])*invMass;

                                                                    pos_ptr[index].x=pos_ptr[index].x+dt_temp*vel_ptr[index].x;
                                                                    pos_ptr[index].y=pos_ptr[index].y+dt_temp*vel_ptr[index].y;
                                                                    pos_ptr[index].z=pos_ptr[index].z+dt_temp*vel_ptr[index].z;

                                                                    force_ptr[index]=make_real4(0.0);});
                        CudaSafeCall(cudaStreamSynchronize(stream));
                    }

                    initialized = true;
                }

                void forwardTime() override {

                    if(not initialized){
                        this->init();
                    }

                    System::log<System::DEBUG1>("[BBK] (%s) Performing integration step %llu",name.c_str(), this->gd->getFundamental()->getCurrentStep());

                    if(stopTransRotSteps > 0 ){
                       if(this->gd->getFundamental()->getCurrentStep()%stopTransRotSteps == 0){
                           this->copyToRef();
                       }
                    }

                    this->updateForce();
                    this->integrationStep();

                    if(stopTransRotSteps > 0){
                        if(this->gd->getFundamental()->getCurrentStep()%stopTransRotSteps == 0){
                            //Stop global movements
                            this->stopTransRot();
                        }
                    }

                    this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
                    this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
                    CudaSafeCall(cudaStreamSynchronize(stream));
                }

                void integrationStep(){

                    int N = this->pg->getNumberParticles();

                    auto mass     = this->pd->getMass(access::location::gpu, access::mode::read);
                    auto friConst = this->pd->getFrictionConstant(access::location::gpu, access::mode::read);

                    auto pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                    auto vel   = this->pd->getVel(access::location::gpu, access::mode::readwrite);
                    auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                    auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                    real* mass_ptr = mass.raw();
                    real* friConst_ptr = friConst.raw();

                    real4* pos_ptr   = pos.raw();
                    real3* vel_ptr   = vel.raw();
                    real4* force_ptr = force.raw();

                    uint step_temp = this->gd->getFundamental()->getCurrentStep();
                    uint seed_temp = this->gd->getSystem()->getSeed();
                    real dt_temp   = this->dt;
                    real kBT_temp  = this->kBT;
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

                    CudaCheckError();

                }
        };
    }

}}}

#endif

