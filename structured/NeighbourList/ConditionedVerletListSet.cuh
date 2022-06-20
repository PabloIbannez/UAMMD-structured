#ifndef CONDITIONEDVERLETLISTSET_CUH
#define CONDITIONEDVERLETLISTSET_CUH

#include"./ConditionedListSetBase.cuh"

#include<limits>

namespace uammd{
namespace structured{

  namespace ConditionedVerletListSet_ns{
    
    template<class PosIterator>
    __global__ void checkMaximumDrift(PosIterator currentPos,
				                      real4 *storedPos,
				                      real maxDistAllowed,
				                      uint* errorFlag,
				                      Box box,
				                      int numberParticles){
      int id = blockIdx.x*blockDim.x + threadIdx.x;
      if(id>=numberParticles) return;
      
      real3 currentPos_i  = make_real3(currentPos[id]);
      real3 previousPos_i = make_real3(storedPos[id]);
      
      real3 rij = box.apply_pbc(currentPos_i - previousPos_i);

      if(dot(rij, rij)>=(maxDistAllowed*maxDistAllowed)){atomicAdd(errorFlag, 1u);}
    }
  
  }

  template<class condition>
  class ConditionedVerletListSet: public ParameterUpdatable{
      
      private:
          
          shared_ptr<System> sys;
          shared_ptr<ParticleData> pd;
          shared_ptr<ParticleGroup> pg;
          
          thrust::device_vector<real4> storedPos;
          
          uint*                        errorFlagsCPU;
          thrust::device_vector<uint>  errorFlagsGPU;
          
          ConditionedListSetBase<condition> nl;
            
          connection posWriteConnection, reorderConnection;
          
          real verletRadiusMultiplier;
          int stepsSinceLastUpdate = 0;
          
          bool forceNextRebuild=true;

          Box box;
          
          real cutOff;
          real cutOffVerlet;
                
          real thresholdDistance;
      
      public:

          struct Parameters{
              real cutOff;
              real cutOffVerlet;
          };

          using NeighbourListData = typename ConditionedListSetBase<condition>::NeighbourListData;

          ConditionedVerletListSet(shared_ptr<ParticleGroup> pg,
                                   shared_ptr<condition>   cond,
                                   Parameters par):pg(pg),
                                                   pd(pg->getParticleData()), 
                                                   sys(pg->getParticleData()->getSystem()), 
                                                   nl(ConditionedListSetBase<condition>(pg,cond)),
                                                   cutOff(par.cutOff),
                                                   cutOffVerlet(par.cutOffVerlet){

              sys->log<System::MESSAGE>("[ConditionedVerletListSet] Created");

              cudaError_t status = cudaMallocHost((void**)&errorFlagsCPU, sizeof(uint));
              if (status != cudaSuccess){
                  sys->log<System::CRITICAL>("[ConditionedVerletListSet] Error allocating pinned host memory");
              }
              
              errorFlagsGPU.resize(1);
              
              reorderConnection  = pd->getReorderSignal()->connect([this](){this->handleReorder();});
            
              thresholdDistance = (cutOffVerlet-cutOff)/2.0;
              
              CudaCheckError();
          }

          ~ConditionedVerletListSet(){
                posWriteConnection.disconnect();
                reorderConnection.disconnect();
                cudaFreeHost(errorFlagsCPU);
          }
          
          NeighbourListData getNeighbourList(std::string conditionName,cudaStream_t st = 0){
              auto   listData = nl.getNeighbourList(conditionName,st);
              return listData;
          }
          
          void update(cudaStream_t st = 0){
              if(needsRebuild(box,st)){
                  //System::log<System::MESSAGE>("[VerletList] Verlet list rebuild. Last one %d steps ago.",stepsSinceLastUpdate);
                  stepsSinceLastUpdate = 0;
                
                  storeCurrentPos(st);
                  rebuildList(box,st);
              }

              stepsSinceLastUpdate++;
          }
          
          int getNumberOfStepsSinceLastUpdate(){
              return stepsSinceLastUpdate-1;
          }

          void updateBox(Box newBox){
              box = newBox;
          }

          real getCutOff(){return cutOff;}
          real getCutOffVerlet(){return cutOffVerlet;}

          void setCutOff(real newCutOff){ 
              this->cutOff = newCutOff;
              thresholdDistance = (cutOffVerlet-cutOff)/2.0;
              forceNextRebuild = true;
          }
          
          void setCutOffVerlet(real newCutOffVerlet){ 
              this->cutOffVerlet = newCutOffVerlet;
              thresholdDistance = (cutOffVerlet-cutOff)/2.0;
              forceNextRebuild = true;
          }

      private:
            
          void handleReorder(){
              sys->log<System::DEBUG1>("[ConditionedVerletListSet] Issuing a list update after a reorder.");
              forceNextRebuild = true;
          }

          void storeCurrentPos(cudaStream_t st = 0){
                auto pos = pd->getPos(access::location::gpu, access::mode::read);
                auto posIter = pg->getPropertyIterator(pos);
                      
                int numberParticles = pg->getNumberParticles();
                
                storedPos.resize(numberParticles);
                int Nthreads=256;
                //int Nthreads=1024;
                int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);
                copyGPU<<<Nthreads,Nblocks,0,st>>>(posIter,
                                                   thrust::raw_pointer_cast(storedPos.data()),
                                                   numberParticles);
                CudaCheckError();
          }

          void rebuildList(Box box,cudaStream_t st){

            auto pos = pd->getPos(access::location::gpu, access::mode::read);
            auto posIter = pg->getPropertyIterator(pos);
                  
            int numberParticles = pg->getNumberParticles();

            nl.update(posIter, numberParticles, box, cutOffVerlet, st);
          }

          bool needsRebuild(Box box,cudaStream_t st = 0){
              
              if(forceNextRebuild){
                  forceNextRebuild = false;
                  return true;
              }
              
              if(isParticleDriftOverThreshold(box,st)){
                  System::log<System::DEBUG2>("[VerletList] Particle drift forced rebuild. Last one %d steps ago.",
                          stepsSinceLastUpdate);
                  return true;
              }

              return false;
          }

          bool isParticleDriftOverThreshold(Box box,cudaStream_t st = 0){
              
              auto pos = pd->getPos(access::location::gpu, access::mode::read);
              auto posIter = pg->getPropertyIterator(pos);
                    
              int numberParticles = pg->getNumberParticles();

              errorFlagsCPU[0] = 0;
              cudaMemcpyAsync(thrust::raw_pointer_cast(errorFlagsGPU.data()), errorFlagsCPU, 
                              sizeof(uint), cudaMemcpyHostToDevice,st);
              
              int blockSize = 256;
              int nblocks = numberParticles/blockSize+1;

              ConditionedVerletListSet_ns::checkMaximumDrift<<<nblocks, blockSize, 0, st>>>(posIter,
                                                                                            thrust::raw_pointer_cast(storedPos.data()),
                                                                                            thresholdDistance,
                                                                                            thrust::raw_pointer_cast(errorFlagsGPU.data()),
                                                                                            box,
                                                                                            numberParticles);
              CudaCheckError();
              cudaMemcpyAsync(errorFlagsCPU,
                              thrust::raw_pointer_cast(errorFlagsGPU.data()),  
                              sizeof(uint), cudaMemcpyDeviceToHost,st);

              cudaStreamSynchronize(st);
              uint errorFlag = errorFlagsCPU[0];
              if(errorFlag>0){
                  System::log<System::DEBUG2>("[VerletList] Found %d particles over threshold.", errorFlag);
              }
              bool   isSomeParticleDisplacementOverThreshold = errorFlag>0;
              return isSomeParticleDisplacementOverThreshold;
          }

  };

}}
#endif
