#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

    namespace HessianMeasure_ns{

      void resetForces(std::shared_ptr<ParticleData> pd, cudaStream_t st){
	auto force = pd->getForce(access::location::gpu, access::mode::write);
	thrust::fill(thrust::cuda::par.on(st), force.begin(), force.end(), real4());
	CudaCheckError();
      }

      // Displaces the particle specified by 'idParticle' by a distance 'dr' in the direction
      // specified by 'dir' (where 0,1,2 corresponds to x,y,z respectively).
      void moveElement(std::shared_ptr<ParticleData> pd, real dr, int idParticle, int dir){
	auto pos = pd->getPos(access::location::cpu, access::mode::readwrite);
	if (dir == 0)
	  pos[idParticle].x += dr;
	else if (dir == 1)
	  pos[idParticle].y += dr;
	else if (dir == 2)
	  pos[idParticle].z += dr;
      }

      // Computes the forces acting on all particles by iterating over a list of interactors.
      thrust::host_vector<real4> computeForces(std::shared_ptr<ParticleData> pd,
					       std::map<std::string,std::shared_ptr<typename uammd::Interactor>> interactors,
					       cudaStream_t st){
	// Reset forces before computing
	resetForces(pd, st);

	// Sum the forces from each interactor
	for(auto interactor: interactors)
	  interactor.second->sum({.force =true, .energy = false, .virial = false}, st);
	CudaCheckError();

	// Fetch the computed forces to host memory
	auto forcePD = pd->getForce(access::location::cpu, access::mode::read);
	thrust::host_vector<real4> force_cpu(forcePD.size());
	thrust::copy(forcePD.begin(), forcePD.end(), force_cpu.begin());
	return force_cpu;
      }

      void computeHessian_i_Numerically(std::shared_ptr<ParticleData> pd,
					std::map<std::string,std::shared_ptr<typename uammd::Interactor>> interactors,
					int index_particle, cudaStream_t st){

	int nParticles = pd->getNumParticles();
	real dr        = sqrt(std::numeric_limits<real>::epsilon());
	auto hessian_i = pd->getHessian(access::location::cpu, access::mode::readwrite);

	// Reference forces without any displacement
	auto f_ref = computeForces(pd, interactors, st);

	//x direction
	moveElement(pd, dr, index_particle, 0);
	auto forces_i_moved_x = computeForces(pd, interactors, st);
	forj(0, nParticles){
	  real3 df = make_real3(forces_i_moved_x[j] - f_ref[j])/dr;
	  hessian_i[j].xx = -df.x;
	  hessian_i[j].xy = -df.y;
	  hessian_i[j].xz = -df.z;
	}
	moveElement(pd, -dr, index_particle, 0);

	//y direction
	moveElement(pd,  dr, index_particle, 1);
	auto forces_i_moved_y = computeForces(pd, interactors, st);
	forj(0, nParticles){
	  real3 df = make_real3(forces_i_moved_y[j] - f_ref[j])/dr;
	  hessian_i[j].yx = -df.x;
	  hessian_i[j].yy = -df.y;
	  hessian_i[j].yz = -df.z;
	}
	moveElement(pd,  -dr, index_particle, 1);

	//z direction
	moveElement(pd,  dr, index_particle, 2);
	auto forces_i_moved_z = computeForces(pd, interactors, st);
	forj(0, nParticles){
	  real3 df = make_real3(forces_i_moved_z[j] - f_ref[j])/dr;
	  hessian_i[j].zx = -df.x;
	  hessian_i[j].zy = -df.y;
	  hessian_i[j].zz = -df.z;
	}
	moveElement(pd,  -dr, index_particle, 2);
      }
    }

  class HessianMeasure: public SimulationStepBase_EnergyForceTorqueHessian{

    std::string    outputFilePath;
    std::ofstream  outputFile;
    std::string mode;
    bool firstMeasure;
    int outputPrecision;

  public:

    HessianMeasure(std::shared_ptr<ParticleGroup>             pg,
		   std::shared_ptr<IntegratorManager> integrator,
		   std::shared_ptr<ForceField>                ff,
		   DataEntry& data,
		   std::string name):SimulationStepBase_EnergyForceTorqueHessian(pg,
							                         integrator,ff,
							                         data,name){

      //Read parameters
      outputFilePath  = data.getParameter<std::string>("outputFilePath");
      mode            = data.getParameter<std::string>("mode");
      outputPrecision = data.getParameter<int>("outputPrecision",10);
      firstMeasure    = true;
    }


    void init(cudaStream_t st) override{

        bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

        if (isFileEmpty){
            outputFile<<"#"<< std::left << std::fixed
                      << " "  << "Id_i"<< " "  << "Id_j"
                      << " " << "Hxx" << " " << "Hxy" << " "  << "Hxz"
                      << " " << "Hyx" << " " << "Hyy" << " "  << "Hyz"
                      << " " << "Hzx" << " " << "Hzy" << " "  << "Hzz"
                      <<"\n";

        } else {
            firstMeasure = false;
        }
    }

    void applyStep(ullint step, cudaStream_t st) override{

      //Set hessian to zero
      auto groupIndex     = pg->getIndexIterator(access::location::gpu);
      int numberParticles = pg->getNumberParticles();

      auto id = this->pd->getId(access::location::cpu, access::mode::read);

      if (!firstMeasure){
	outputFile<<"#\n";
      } else {
	firstMeasure = false;
      }

      fori(0,numberParticles){
	int index_i    = groupIndex[i];
        int id_i       = id[index_i];
	{
	  auto hessian_i = this->pd->getHessian(access::location::gpu, access::mode::write);
	  auto selId     = this->pd->getSelectedId(access::location::gpu, access::mode::write);
	  thrust::fill(thrust::cuda::par.on(st),
		       hessian_i.begin(),
		       hessian_i.end(),
		       tensor3());


	  thrust::fill(thrust::cuda::par.on(st),
		       selId.begin(),
		       selId.end(),
		       id_i);
	}
	//Sum hessian
	if (mode == "Analytical"){
	  for(auto& interactor : this->topology->getInteractors()){
	    //Create computable
	    uammd::Interactor::Computables compTmp;
	    compTmp.hessian = true;
	    interactor.second->sum(compTmp,st);
	  }
	} else if (mode == "Numerical"){
	  HessianMeasure_ns::computeHessian_i_Numerically(pd, this->topology->getInteractors(),
							  index_i, st);
	} else {
	  System::log<System::CRITICAL>("[HessianMeasure] Unknown mode: ", mode);
        }

	auto hessian_i = this->pd->getHessian(access::location::cpu, access::mode::read);
	forj(0,numberParticles){
	  int index_j = groupIndex[j];
	  int id_j    = id[index_j];
	  tensor3 Hij = hessian_i[index_j];

	  outputFile << std::left << std::fixed << std::setprecision(outputPrecision)
		     << " " << id_i   << " " << id_j
		     << " " << Hij.xx << " " << Hij.xy << " " << Hij.xz
		     << " " << Hij.yx << " " << Hij.yy << " " << Hij.yz
		     << " " << Hij.xz << " " << Hij.yz << " " << Hij.zz;
	  outputFile<<"\n";
	}
      }
    }
  };

}}}}

REGISTER_SIMULATION_STEP(
    MechanicalMeasure,HessianMeasure,
    uammd::structured::SimulationStep::SimulationMeasures::HessianMeasure
)
