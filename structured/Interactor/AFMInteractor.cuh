#ifndef __AFM_INTERACTOR__
#define __AFM_INTERACTOR__

#include <thrust/iterator/discard_iterator.h>

namespace uammd{
namespace structured{

    namespace Interactor{
    namespace SphericalTip_ns{

        void __global__ computeForceKernel(real4* pos,
                                           real4* force,
                                           real4* tipForce,
                                           int numberParticles,
			                               ParticleGroup::IndexIterator groupSampleIterator,
			                               ParticleGroup::IndexIterator groupTipIterator,
                                           real epsilonTip,
                                           real sigmaTip,
                                           real A,real B,
                                           real epsilonTipSurf,
                                           real sigmaTipSurf,
                                           real ATipSurf,real BTipSurf,
                                           real Rtip,
                                           real3 chipPos,
                                           real  surfPos,
                                           real Kxy,
                                           real K){
            
            int i = blockIdx.x*blockDim.x + threadIdx.x;
            
            const real3 tPos = make_real3(pos[groupTipIterator[0]]);
            
            if(i < numberParticles){

                const real3 p    = make_real3(pos[groupSampleIterator[i]]);

                const real3 dr = tPos-p;

                const real r  = sqrt(dot(dr,dr));
                      real r2 = r-Rtip;
                           r2 = r2*r2;

                const real sinvr2  = sigmaTip*sigmaTip/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -epsilonTip*real(6.0)*(real(2.0)*A*sinvr12+B*sinvr6)/abs(r-Rtip);

                real3 f = fmod*(dr/r);

                //Force on the particles due to the tip
                force[groupSampleIterator[i]] += make_real4(f,0);
                
                //Force on the tip due to the particles and the umbrella-like potential
                tipForce[i] = make_real4(0,0,-f.z,0); //Note NOT +=
            }
            
            //Compute surface and umbrella force in thread 0
            if(i==0){

                real fmodSurf = real(0.0);
                    
                const real dz  = abs(surfPos-tPos.z);

                if((ATipSurf == real(0.0)) and (BTipSurf == real(0.0))){
                
                } else {
                    
                    real dz2 = dz-Rtip;
                         dz2 = dz2*dz2;
                    
                    const real sinvdz2  = sigmaTipSurf*sigmaTipSurf/dz2;
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;
                    
                    fmodSurf = -epsilonTipSurf*real(6.0)*(real(2.0)*ATipSurf*sinvdz12+BTipSurf*sinvdz6)/abs(dz-Rtip); 
                }
                
                if(numberParticles == 0){
                    tipForce[0] = make_real4(0);
                }

                
                tipForce[0] += make_real4(Kxy*(chipPos.x-tPos.x),
                                          Kxy*(chipPos.y-tPos.y),
                                          K*(chipPos.z-tPos.z)+fmodSurf*(surfPos-tPos.z)/dz,
                                          0);
            }
            
        }
    }
    
    class SphericalTip: public Interactor{

        private:

            int N;

            std::shared_ptr<ParticleGroup> pgTip;
                
            thrust::device_vector<real4> tipForceTemp;
                
            void*  cubReductionTempStorage    ;
            size_t cubReductionTempStorageSize;
            
            double3 chipPos;
            double  surfPos;
            
        public:

            struct Parameters{
                
                real epsilonTip;
                real sigmaTip;
                real A;
                real B;
                
                real epsilonTipSurf;
                real sigmaTipSurf;
                real ATipSurf;
                real BTipSurf;
                
                real Rtip;

                real Kxy;
                real K;
            };

        private:
            
            real epsilonTip;
            real sigmaTip;
            real A;
            real B;
            
            real epsilonTipSurf;
            real sigmaTipSurf;
            real ATipSurf;
            real BTipSurf;
            
            real Rtip;
            
            real Kxy;
            real K;

        public:
            
            Parameters inputFileToParam(InputFile& in){

                Parameters param;
                
                in.getOption("epsilonTip",InputFile::Required)
                              >>param.epsilonTip;
                in.getOption("sigmaTip",InputFile::Required)
                              >>param.sigmaTip;

                in.getOption("Atip",InputFile::Required)
                              >>param.A;
                in.getOption("Btip",InputFile::Required)
                              >>param.B;
                
                in.getOption("epsilonTipSurf",InputFile::Required)
                              >>param.epsilonTipSurf;
                in.getOption("sigmaTipSurf",InputFile::Required)
                              >>param.sigmaTipSurf;
                in.getOption("ATipSurf",InputFile::Required)
                              >>param.ATipSurf;
                in.getOption("BTipSurf",InputFile::Required)
                              >>param.BTipSurf;
                
                in.getOption("Rtip",InputFile::Required)
                              >>param.Rtip;
                
                in.getOption("Ktip",InputFile::Required)
                              >>param.K;
                in.getOption("Kxytip",InputFile::Required)
                              >>param.Kxy;


                return param;

            }

            SphericalTip(std::shared_ptr<System>       sys,
                         std::shared_ptr<ParticleData>  pd,
                         std::shared_ptr<ParticleGroup> pgSample,
                         std::shared_ptr<ParticleGroup> pgTip,
                         InputFile& in):SphericalTip(sys,pd,pgSample,pgTip,inputFileToParam(in)){}

            SphericalTip(std::shared_ptr<System>       sys,
                         std::shared_ptr<ParticleData>  pd,
                         std::shared_ptr<ParticleGroup> pgSample,
                         std::shared_ptr<ParticleGroup> pgTip,
                         Parameters par):Interactor(pd,pgSample,sys,std::string("SphericalTip")),
                                         pgTip(pgTip),
                                         epsilonTip(par.epsilonTip),
                                         sigmaTip(par.sigmaTip),
                                         A(par.A),
                                         B(par.B),
                                         epsilonTipSurf(par.epsilonTipSurf),
                                         sigmaTipSurf(par.sigmaTipSurf),
                                         ATipSurf(par.ATipSurf),
                                         BTipSurf(par.BTipSurf),
                                         Rtip(par.Rtip),
                                         Kxy(par.Kxy),
                                         K(par.K)
            {
                
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter epsilonTip %f",epsilonTip);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter sigmaTip %f",sigmaTip);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter Atip %f",A);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter Btip %f",B);
                
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter epsilonTipSurf %f",epsilonTipSurf);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter sigmaTipSurf %f",sigmaTipSurf);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter ATipSurf %f",ATipSurf);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter BTipSurf %f",BTipSurf);
                
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter Rtip %f",Rtip);

                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter Kxytip %f",Kxy);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter Ktip %f",K);
                
                this->N=this->pg->getNumberParticles();
                if(this->N!=0){
                    tipForceTemp.resize(this->N);
                } else {
                    tipForceTemp.resize(1);
                }

                cubReductionTempStorage     = NULL;
                cubReductionTempStorageSize = 0;
                
                auto tipForceIter = pgTip->getPropertyIterator(pd->getForce(access::location::gpu, access::mode::readwrite));
                
                if(this->N==0){
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::raw_pointer_cast(tipForceTemp.data()), 
                                           tipForceIter,
                                           1);
                } else {
                    cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                           cubReductionTempStorageSize,
                                           thrust::raw_pointer_cast(tipForceTemp.data()), 
                                           tipForceIter,
                                           this->N);
                }

                cudaMalloc(&cubReductionTempStorage, cubReductionTempStorageSize);

                cudaDeviceSynchronize();
            }
            
            ~SphericalTip(){
                cudaFree(cubReductionTempStorage);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                
                int blocksize = 256;
                
                int Nthreads;
                int Nblocks;
                    
                if(this->N!=0){
                    Nthreads = (blocksize<this->N)?blocksize:this->N;
                    Nblocks=this->N/Nthreads + ((this->N%Nthreads)?1:0);
                } else {
                    Nthreads = 1;
                    Nblocks  = 1;
                }
                    
                if(comp.force == true){
                
                    auto pos   = pd->getPos(access::location::gpu, access::mode::readwrite).raw();
                    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();

                    auto tipForceTemp_ptr = thrust::raw_pointer_cast(tipForceTemp.data());

                    auto gSampleIter = pg->getIndexIterator(access::location::gpu);
                    auto gTipIter    = pgTip->getIndexIterator(access::location::gpu);

                    SphericalTip_ns::computeForceKernel<<<Nblocks, Nthreads, 0, st>>>(pos,
                                                                                      force,tipForceTemp_ptr,
                                                                                      this->N,
                                                                                      gSampleIter,
                                                                                      gTipIter,
                                                                                      epsilonTip,
                                                                                      sigmaTip,
                                                                                      A,B,
                                                                                      epsilonTipSurf,
                                                                                      sigmaTipSurf,
                                                                                      ATipSurf,BTipSurf,
                                                                                      Rtip,
                                                                                      {real(chipPos.x),real(chipPos.y),real(chipPos.z)},
                                                                                      real(surfPos),
                                                                                      Kxy,
                                                                                      K);  
                    //CudaSafeCall(cudaStreamSynchronize(st));
                    
                    auto tipForceIter = pgTip->getPropertyIterator(pd->getForce(access::location::gpu, access::mode::readwrite));
                    
                    if(this->N!=0){
                        cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                               cubReductionTempStorageSize,
                                               tipForceTemp_ptr, 
                                               tipForceIter,
                                               this->N,st);
                    } else {
                        cub::DeviceReduce::Sum(cubReductionTempStorage    ,
                                               cubReductionTempStorageSize,
                                               tipForceTemp_ptr, 
                                               tipForceIter,
                                               1,st);
                    }

                }

                if(comp.energy == true){

                }
                
                if(comp.virial == true){

                }


            }
            
            void    setChipPosition(double3 newChipPos){chipPos=newChipPos;}
            double3 getChipPosition(){return chipPos;}
            
            void    setChipHeight(double newChipHeight){
                                    this->setChipPosition({chipPos.x,chipPos.y,newChipHeight});
            }
            double  getChipHeight(){return chipPos.z;}
            
            void   setSurfacePosition(double newSurfPos){surfPos=newSurfPos;}
            double getSurfacePosition(){return surfPos;}

            real3 getTipPosition(){
                    
                auto pos   = pd->getPos(access::location::cpu, access::mode::read);
                auto gTipIter = pgTip->getIndexIterator(access::location::cpu);
                return make_real3(pos[gTipIter[0]]);
            }
            
            real getTipDeflection(){
                return (this->getTipPosition().z-this->getChipPosition().z);
            }

            real getTipDeflectionForce(){
                return K*this->getTipDeflection();
            }

            template<class UNITS>
            void applyUnits(){

                this->Kxy = this->Kxy*UNITS::TO_INTERNAL_FORCE; 
                this->K   = this->K*UNITS::TO_INTERNAL_FORCE;
                
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter Kxytip (after units): %f",Kxy);
                this->sys->log<System::MESSAGE>("[SphericalTip] "
                                                 "Parameter Ktip (after units): %f",K);

            }
    };

}}}

#endif
