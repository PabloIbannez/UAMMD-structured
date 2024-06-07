#ifndef __PAIRWISEFORCE_MEASURE__
#define __PAIRWISEFORCE_MEASURE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

  class PairwiseForceMeasure: public SimulationStepBase_EnergyForceTorqueHessian{

    std::string    outputFilePath;
    std::ofstream  outputFile;
    std::string    mode;
    bool firstMeasure;
    std::vector<std::string> availableModes = {"Total_force", "Pairwise_force"};
  public:

    PairwiseForceMeasure(std::shared_ptr<ParticleGroup>             pg,
			 std::shared_ptr<IntegratorManager> integrator,
			 std::shared_ptr<ForceField>                ff,
			 DataEntry& data,
			 std::string name):SimulationStepBase_EnergyForceTorqueHessian(pg,
										       integrator,ff,
										       data,name){

      //Read parameters
      outputFilePath = data.getParameter<std::string>("outputFilePath");
      mode           = data.getParameter<std::string>("mode","Pairwise_force");
      firstMeasure   = true;
      bool validMode = false;
      fori(0,availableModes.size()){
	if (mode == availableModes[i]){
	  validMode = true;
	  break;
	}
      }
      if (!validMode){
	System::log<System::CRITICAL>("[PairwiseForceMeasure] Unknown mode: %s", mode);
      }
    }

    void init(cudaStream_t st) override{
      bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

      if (isFileEmpty){
        outputFile<<"#"<< std::left << std::fixed;

          if (mode == "Pairwise_force"){
        	  outputFile << std::setw(7) << "From" << std::setw(8)  << "To";
          } else {
            outputFile << std::setw(7) << "Id";
          }

        outputFile << std::setw(20)  << "Fx" << std::setw(20)  << "Fy" << std::setw(20)  << "Fz";
        outputFile <<"\n";
      } else {
        firstMeasure = false;
      }
    }

    void applyStep(ullint step, cudaStream_t st) override{

      //Set hessian to zero
      auto groupIndex     = pg->getIndexIterator(access::location::gpu);
      int numberParticles = pg->getNumberParticles();
      auto id             = this->pd->getId(access::location::cpu, access::mode::read);

      if (!firstMeasure){
	outputFile<<"#\n";
      } else {
	firstMeasure = false;
      }

      if (mode == "Pairwise_force"){
	fori(0,numberParticles){
	  int index_i    = groupIndex[i];
	  int id_i       = id[index_i];
	  {
	    auto pairforce_i = this->pd->getPairwiseForce(access::location::gpu, access::mode::write);
	    //selId is the id of the particle that recieves the forces from the other particles
	    auto selId       = this->pd->getSelectedId(access::location::gpu, access::mode::write);
	    thrust::fill(thrust::cuda::par.on(st),
			 pairforce_i.begin(),
			 pairforce_i.end(),
			 real4());


	    thrust::fill(thrust::cuda::par.on(st),
			 selId.begin(),
			 selId.end(),
			 id_i);
	  }
	  //Compute pairwise forces
	  for(auto& interactor : this->topology->getInteractors()){
	    //Create computable
	    uammd::Interactor::Computables compTmp;
	    compTmp.pairwiseForce = true;
	    interactor.second->sum(compTmp,st);
	  }

	  auto pairforces_i = this->pd->getPairwiseForce(access::location::cpu, access::mode::read);
	  forj(0,numberParticles){
	    int index_j = groupIndex[j];
	    int id_j    = id[index_j];
	    real4 f_ij  = pairforces_i[index_j];
	    if (f_ij.x != 0.0 or f_ij.y != 0.0 or f_ij.z != 0.0){
	      outputFile << std::left << std::fixed
			 << std::setprecision(10)
			 << std::setw(8)  << id_j
			 << std::setw(8)  << id_i
			 << std::setw(20) << f_ij.x
			 << std::setw(20) << f_ij.y
			 << std::setw(20) << f_ij.z <<"\n";
	    }
	  }
	}
      } else if (mode == "Total_force"){
	{
	  auto forces = this->pd->getPairwiseForce(access::location::gpu, access::mode::write);
	  thrust::fill(thrust::cuda::par.on(st),
		       forces.begin(),
		       forces.end(),
		       real4());
	}
	//Compute forces
	for(auto& interactor : this->topology->getInteractors()){
	  //Create computable
	  uammd::Interactor::Computables compTmp;
	  compTmp.force = true;
	  interactor.second->sum(compTmp,st);
	}

	auto forces = this->pd->getForce(access::location::cpu, access::mode::read);
	forj(0,numberParticles){
	  int index_j = groupIndex[j];
	  int id_j    = id[index_j];
	  real4 f_j = forces[index_j];
	  if (f_j.x != 0.0 or f_j.y != 0.0 or f_j.z != 0.0){
	    outputFile << std::left << std::fixed
		       << std::setprecision(10)
		       << std::setw(8)  << id_j
		       << std::setw(20) << f_j.x
		       << std::setw(20) << f_j.y
		       << std::setw(20) << f_j.z <<"\n";
	  }
	}
      }
    }
  };

}}}}

#endif
