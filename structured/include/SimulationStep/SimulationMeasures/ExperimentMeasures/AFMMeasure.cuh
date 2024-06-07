#ifndef __AFM_MEASURE__
#define __AFM_MEASURE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class AFMMeasure: public SimulationStepBase_EnergyForceTorque{

  int batchId;

  std::shared_ptr<Interactor::AFMInteractor<Potentials::AFM::SphericalTip>> afm;

  std::shared_ptr<ParticleGroup> tip;
  std::shared_ptr<ParticleGroup> sample;
  //

  std::string   outputFilePath;
  std::ofstream outputFile;

  bool verbose;

  std::string tipType;

public:

  AFMMeasure(std::shared_ptr<ParticleGroup>              pg,
			       std::shared_ptr<IntegratorManager>  integrator,
			       std::shared_ptr<ForceField> ff,
			       DataEntry& data,
			       std::string name):SimulationStepBase_EnergyForceTorque(pg,integrator,ff,data,name){

    outputFilePath          = data.getParameter<std::string>("outputFilePath");

    verbose                 = data.getParameter<bool>("verbose",false);

    //

    tipType = data.getParameter<std::string>("tipType","TIP");

    System::log<System::MESSAGE>("[AFMMeasure] Output file path: %s" , outputFilePath.c_str());
    System::log<System::MESSAGE>("[AFMMeasure] Verbose: %s" , verbose ? "true" : "false");

  }

  ~AFMMeasure(){
    outputFile.close();
  }

  void init(cudaStream_t st) override{

    System::log<System::MESSAGE>("[AFMMeasure] Initializing...");

    int batchNumber = GroupUtils::BatchGroupNumber(this->pg);
		if(batchNumber > 1){
			System::log<System::CRITICAL>("[AFMMMeasure] Only one batch is allowed");
		} else {
      batchId = GroupUtils::getBatchesInGroup(this->pg)[0];
    }

    //Open output file

    bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

    if(isFileEmpty){
      if(verbose){
        outputFile << std::setw(24) << "#Step ";
        outputFile << std::setw(24) << "ChipPos ";
        outputFile << std::setw(24) << "TipPos ";
        outputFile << std::setw(24) << "tipDeflection ";
        outputFile << std::setw(24) << "tipForce ";
        outputFile << std::setw(24) << "sampleForce ";
        outputFile << std::setw(24) << "tipDeflectionForce " << std::endl;
      } else {
        outputFile << std::setw(24) << "#Step ";
        outputFile << std::setw(24) << "X ";
        outputFile << std::setw(24) << "F " << std::endl;
      }
    }

    //Get AFM interactor
		for(auto& interactor : this->topology->getInteractors()){
      std::string interactorName = interactor.first;

      InputEntryManager::entryInfo interactorInfo = this->topology->getInteractorInfo(interactorName);

      std::string type    = interactorInfo.entryType;
      std::string subType = interactorInfo.entrySubType;

      if(type == "AFM"){
        //Check if afm interactor is already set
        if(afm){
          System::log<System::CRITICAL>("[AFMMeasure] More than one AFM interactor found!");
        } else {
          if(subType == "SphericalTip"){
            afm = std::dynamic_pointer_cast<Interactor::AFMInteractor<Potentials::AFM::SphericalTip>>(interactor.second);
          } else {
            System::log<System::CRITICAL>("[AFMMeasure] Unknown AFM interactor subtype!");
          }
        }
      }
		}

    //Create tip and sample groups
    auto types = this->gd->getTypes();
    std::vector<int> tipTypeId;
    tipTypeId.push_back(types->getTypeId(tipType));

    GroupUtils::selectors::types tipSelector(tipTypeId);
    tip    = GroupUtils::createSubGroup(this->pg,tipSelector,pg->getName() + "_tip");

    GroupUtils::selectors::notTypes sampleSelector(tipTypeId);
    sample = GroupUtils::createSubGroup(this->pg,sampleSelector,pg->getName() + "_sample");

    //Check tip has only one particle
    if(tip->getNumberParticles() != 1){
      System::log<System::CRITICAL>("[AFMMeasure] Tip group must have only one particle!");
    }

    System::log<System::MESSAGE>("[AFMMeasure] Initializing end");
  }

  void applyStep(ullint step, cudaStream_t st) override {

    //Calling SimulationStepBase_EnergyForceTorque::setZero,
    //the energy, force and torque are set to zero for the whole group
    this->setZero(st);

    //////////////////////////////////////////////////

    auto afmParams = afm->getAFMParameters()[batchId];

    real dt = this->gd->getFundamental()->getTimeStep();

    real3  chipPos  = afmParams.startChipPosition;
    ullint indStart = afmParams.indentationStartStep;
    ullint bckwStep = afmParams.indentationBackwardStep;

    if(step >= indStart and step < bckwStep){
        real time = dt*(step - indStart);
        chipPos.z    += afmParams.tipVelocity*time;
    }
    if(step >= bckwStep){
        real time          = dt*(step - bckwStep);
        real timeUntilBckw = dt*(bckwStep - indStart);
        chipPos.z   += afmParams.tipVelocity*timeUntilBckw - afmParams.tipVelocity*time;
    }

    real3 tipPos;
    {
      auto groupIndex = tip->getIndexIterator(access::location::cpu);
      auto ids = this->pd->getId(uammd::access::location::cpu,uammd::access::mode::read);
      tipPos = make_real3(this->pd->getPos(uammd::access::location::cpu,uammd::access::mode::read)[groupIndex[0]]);
    }

    real tipDeflection = tipPos.z - chipPos.z;

    //////////////////////////////////////////////////
    // Update forces

    {
      uammd::Interactor::Computables comp;
      comp.force  = true;

      this->ff->sum(comp,st);
    }

    // Sum tip force
    real3 tipForce    = Measures::totalForce(tip,st);
    // Sum sample force
    real3 sampleForce = Measures::totalForce(sample,st);

    //////////////////////////////////////////////////
    // Compute tip deflection force

    real K = afmParams.K;
    real tipDeflectionForce = K*tipDeflection;

    //////////////////////////////////////////////////
    // Compute X

    real tipRadius = this->gd->getTypes()->getTypeData(tipType,"radius");
    real X = tipPos.z - tipRadius;

    //////////////////////////////////////////////////
    // Write output

    if(verbose){
      outputFile << std::setw(24) << step;
      outputFile << std::setw(24) << chipPos.z;
      outputFile << std::setw(24) << tipPos.z;
      outputFile << std::setw(24) << tipDeflection;
      outputFile << std::setw(24) << tipForce.z;
      outputFile << std::setw(24) << sampleForce.z;
      outputFile << std::setw(24) << tipDeflectionForce << std::endl;
    } else {
      outputFile << std::setw(24) << step;
      outputFile << std::setw(24) << X;
      outputFile << std::setw(24) << tipDeflectionForce << std::endl;
    }
  }
};

}}}}

#endif
