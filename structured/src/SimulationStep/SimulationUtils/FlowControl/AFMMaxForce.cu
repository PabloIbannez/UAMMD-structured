#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

//TODO: Fix it

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationUtils{

//class AFMMaxForce: public SimulationStepBase{
//
//  std::shared_ptr<Interactor::AFMInteractor<Potentials::AFM::SphericalTip>> afm;
//  std::shared_ptr<ParticleGroup> tip;
//
//  int batchNumber;
//
//  std::map<int,int> batch2tipId;
//  std::vector<bool> isTipMaxForce;
//
//  //
//
//  real maxForce;
//  std::string tipType;
//
//public:
//
//  AFMMaxForce(std::shared_ptr<ParticleGroup>              pg,
//			        std::shared_ptr<IntegratorManager>  integrator,
//			        std::shared_ptr<ForceField> ff,
//			        DataEntry& data,
//			        std::string name):SimulationStepBase(pg,integrator,ff,data,name){
//
//    //
//
//    maxForce = data.getParameter<real>("maxForce");
//
//    tipType = data.getParameter<std::string>("tipType","TIP");
//
//    System::log<System::MESSAGE>("[AFMMaxForce] Max force: %f" , maxForce);
//    System::log<System::MESSAGE>("[AFMMaxForce] Tip type: %s" , tipType.c_str());
//
//  }
//
//  void init(cudaStream_t st) override{
//
//    System::log<System::MESSAGE>("[AFMMaxForce] Initializing...");
//
//    //Get AFM interactor
//		for(auto& interactor : this->topology->getInteractors()){
//      std::string interactorName = interactor.first;
//
//      InputEntryManager::entryInfo interactorInfo = this->topology->getInteractorInfo(interactorName);
//
//      std::string type    = interactorInfo.entryType;
//      std::string subType = interactorInfo.entrySubType;
//
//      if(type == "AFM"){
//        //Check if afm interactor is already set
//        if(afm){
//          System::log<System::CRITICAL>("[AFMMaxForce] More than one AFM interactor found!");
//        } else {
//          if(subType == "SphericalTip"){
//            afm = std::dynamic_pointer_cast<Interactor::AFMInteractor<Potentials::AFM::SphericalTip>>(interactor.second);
//          } else {
//            System::log<System::CRITICAL>("[AFMMaxForce] Unknown AFM interactor subtype!");
//          }
//        }
//      }
//		}
//
//    //Get the number of batches
//    batchNumber = GroupUtils::BatchNumber(this->pd);
//
//    //Create tip
//    auto types = this->gd->getTypes();
//    std::vector<int> tipTypeId;
//    tipTypeId.push_back(types->getTypeId(tipType));
//
//    GroupUtils::selectors::types tipSelector(tipTypeId);
//    tip    = GroupUtils::createSubGroup(this->pg,tipSelector,pg->getName() + "_tip");
//
//    //Check the number of batches matches the number of tip particles
//    if(batchNumber != tip->getNumberParticles()){
//      System::log<System::CRITICAL>("[AFMMaxForce] The number of batches (%d) does not match the number of tip particles (%d)!",
//                                    batchNumber,tip->getNumberParticles());
//    }
//
//    //Create map between batch and tip particle id
//
//    auto ids = this->pd->getId(uammd::access::location::cpu,uammd::access::mode::read);
//    auto groupIndex = tip->getIndexIterator(access::location::cpu);
//    for(int i = 0; i < batchNumber; i++){
//      int tipId = ids[groupIndex[i]];
//      batch2tipId[i] = tipId;
//    }
//
//    //Initialize isTipMaxForce
//    isTipMaxForce.resize(batchNumber,false);
//
//    System::log<System::MESSAGE>("[AFMMaxForce] Initializing end");
//  }
//
//  void applyStep(ullint step, cudaStream_t st) override {
//
//    double time = this->gd->getFundamental()->getSimulationTime();
//    auto id2index = pd->getIdOrderedIndices(access::location::cpu);
//
//    for(int batchId = 0; batchId < batchNumber; batchId++){
//
//      int tipIndex = id2index[batch2tipId[batchId]];
//
//      auto afmParams = afm->getAFMParameters()[batchId];
//
//      real3 chipPos = afmParams.startChipPosition;
//      chipPos.z    += afmParams.tipVelocity*time;
//
//      real3 tipPos;
//      tipPos = make_real3(this->pd->getPos(uammd::access::location::cpu,uammd::access::mode::read)[tipIndex]);
//
//      real tipDeflection      = tipPos.z - chipPos.z;
//
//      real K = afmParams.K;
//      real tipDeflectionForce = K*tipDeflection;
//
//      if(tipDeflectionForce >= maxForce){
//        isTipMaxForce[batchId] = true;
//      }
//    }
//
//    //Check if all tips have reached the max force
//    bool allTipsMaxForce = true;
//    for(int batchId = 0; batchId < batchNumber; batchId++){
//      if(!isTipMaxForce[batchId]){
//        allTipsMaxForce = false;
//        break;
//      }
//    }
//
//    if(allTipsMaxForce){
//      System::log<System::MESSAGE>("[AFMMaxForce] All tips have reached the max force!");
//      this->sys->setState(ExtendedSystem::SIMULATION_STATE::STOPPED);
//    }
//
//  }
//};

}}}}

//REGISTER_SIMULATION_STEP(
//    FlowControl,AfmMaxForce,
//    uammd::structured::SimulationStep::SimulationUtils::AFMMaxForce
//)
