#ifndef __FORCEBETWEENSETS_MEASURE__
#define __FORCEBETWEENSETS_MEASURE__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

  class ForceBetweenSetsMeasure: public SimulationStepBase{

    std::string    outputFilePath;
    std::ofstream  outputFile;
    std::vector<std::string> setNames;
    std::vector<std::vector<int>> particlesInSets;
    bool firstMeasure;

  public:

    ForceBetweenSetsMeasure(std::shared_ptr<ParticleGroup>             pg,
			    std::shared_ptr<IntegratorManager> integrator,
			    std::shared_ptr<ForceField>                ff,
			    DataEntry& data,
			    std::string name):SimulationStepBase(pg,
								 integrator,ff,
								 data,name){

      //Read parameters
      outputFilePath = data.getParameter<std::string>("outputFilePath");
      auto setsData  = data.getDataMap();
      firstMeasure   = true;
      for (auto setData: setsData){
	std::string name = setData.at("name");
	std::vector<int> particlesInSet = setData.at("id_list");
	setNames.push_back(name);
	particlesInSets.push_back(particlesInSet);
      }

      //Check all pariticles in sets belong to the current group
    }



    void init(cudaStream_t st) override{
            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);
            if (isFileEmpty){
                    outputFile<<"#"<< std::left << std::fixed
                            << std::setw(14)  << "From"
                            << std::setw(15)  << "To"
                            << std::setw(20)  << "Fx"
                            << std::setw(20)  << "Fy"
                            << std::setw(20)  << "Fz"<<"\n";
            } else {
                    firstMeasure = false;
            }
    }

    void computePairForcesOverParticle(int id_i, const cudaStream_t &st){
      auto id2index = this->pd->getIdOrderedIndices(access::location::cpu);

      {
	auto pairforce_i = this->pd->getPairwiseForce(access::location::gpu, access::mode::write);
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
    }

    void applyStep(ullint step, cudaStream_t st) override{

      if (!firstMeasure){
	outputFile<<"#\n";
      } else {
	firstMeasure = false;
      }

      fori(0, setNames.size()) {
	std::vector<real3> forces_over_set(setNames.size(), real3());
	for (auto& particleId : particlesInSets[i]) {
	  computePairForcesOverParticle(particleId, st);
	  auto pairforce_i = this->pd->getPairwiseForce(access::location::cpu, access::mode::read);

	  forj(0,setNames.size()){
            if (i == j) continue;
            for (auto& targetParticleId : particlesInSets[j]) {
	      forces_over_set[j] += make_real3(pairforce_i[targetParticleId]);
            }
	  }
	}

	forj(0,setNames.size()) {
	  if (j == i) continue;
	  outputFile << std::left << std::fixed << std::setprecision(10)
		     << std::setw(15) << setNames[i]
		     << std::setw(15) << setNames[j]
		     << std::setw(20) << forces_over_set[j].x
		     << std::setw(20) << forces_over_set[j].y
		     << std::setw(20) << forces_over_set[j].z << "\n";
	}
      }
    }
  };

}}}}

#endif
