//Template: Integrator NVT

#ifndef __INTEGRATOR_MAGNETIC_FORCECOUPLINGMETHOD__
#define __INTEGRATOR_MAGNETIC_FORCECOUPLINGMETHOD__

#include"Integrator/BDHI/BDHI_EulerMaruyama.cuh"
#include "Integrator/BDHI/BDHI_FCM.cuh"

namespace uammd{
namespace structured{
namespace NVT{
namespace Magnetic{

  class ForceCouplingMethod : public IntegratorBasicNVT{
    using Kernel = BDHI::FCM_ns::Kernels::Gaussian;
    using KernelTorque = BDHI::FCM_ns::Kernels::GaussianTorque;
    using FCM_super = BDHI::FCM_impl<Kernel, KernelTorque>;
    
  private:
    std::shared_ptr<FCM_super> fcm;
    std::unique_ptr<MagneticBase> magnetic;
    bool firstStep = true;
    
  public:

    ForceCouplingMethod(std::shared_ptr<GlobalData>           gd,
			std::shared_ptr<ParticleGroup>        pg,
			DataEntry& data,
			std::string name):IntegratorBasicNVT(gd,pg,data,name){


      //Set up the magetic integrator
      typename MagneticBase::Parameters magneticParams;
      
      magneticParams.stream    = this->stream;
      magneticParams.dt        = this->dt;
      magneticParams.kBT       = this->kBT;
      magneticParams.damping   = data.getParameter<real>("damping", -1);
      magneticParams.msat      = data.getParameter<real>("msat");
      magneticParams.gyroRatio = data.getParameter<real>("gyroRatio", -1);
      
      magneticParams.magneticIntegrationAlgorithm       = data.getParameter<std::string>("magneticIntegrationAlgorithm","NotSelected");
      bool resizeBox           = data.getParameter<bool>("resizeBox", false);
      
      magnetic = std::make_unique<MagneticBase>(gd,pg,magneticParams,name);

      //Set up the bdhi integrator
      FCM_super::Parameters bdhiParams;
      bdhiParams.temperature = this ->kBT;
      bdhiParams.viscosity = data.getParameter<real>("viscosity");
      bdhiParams.hydrodynamicRadius = data.getParameter<real>("hydrodynamicRadius", -1.0);
      bdhiParams.tolerance          = data.getParameter<real>("tolerance", 1e-3);
      bdhiParams.cells              = make_int3(data.getParameter<real3>("cells",
									 make_real3(-1, -1, -1)));
      bdhiParams.adaptBoxSize       = data.getParameter<bool>("adaptBoxSize", false);
      bdhiParams.box = gd->getEnsemble()->getBox();
      auto grid      = uammd::BDHI::detail::initializeGrid<Kernel>(bdhiParams);
      bdhiParams.box = grid.box;
      bdhiParams.cells = grid.cellDim;
      bdhiParams.kernel = uammd::BDHI::detail::initializeKernel<Kernel>(bdhiParams, grid);

      
      bdhiParams.hydrodynamicRadius = bdhiParams.kernel->fixHydrodynamicRadius(bdhiParams.hydrodynamicRadius,
									       grid.cellSize.x);
      bdhiParams.kernelTorque = uammd::BDHI::detail::initializeKernelTorque<KernelTorque>(bdhiParams,
											  grid);
      if (bdhiParams.adaptBoxSize){
	gd->getEnsemble()->setBox(grid.box);
	
      }
      fcm = std::make_unique<FCM_super>(bdhiParams);
    }
    
    void updatePosition(){
      int numberParticles = pg->getNumberParticles();
      real dt = gd->getFundamental()->getTimeStep();
      auto indexIter = pg->getIndexIterator(access::location::gpu);
      auto pos = pd->getPos(access::location::gpu, access::mode::readwrite).raw();
      auto dir = pd->getDir(access::location::gpu, access::mode::readwrite).raw();      
      auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
      auto torque = pd->getTorque(access::location::gpu, access::mode::readwrite).raw();    
      auto disp = fcm->computeHydrodynamicDisplacements(pos, force, torque,
							numberParticles, kBT,
							rsqrt(dt), stream);
      auto linearVelocities = disp.first;
      auto angularVelocities = disp.second;
      real3* d_linearV = thrust::raw_pointer_cast(linearVelocities.data());
      real3* d_angularV = thrust::raw_pointer_cast(angularVelocities.data());
      int BLOCKSIZE = 128; 
      int nthreads = BLOCKSIZE<numberParticles?BLOCKSIZE:numberParticles;
      int nblocks = numberParticles/nthreads + ((numberParticles%nthreads!=0)?1:0);      
      uammd::BDHI::FCM_ns::integrateEulerMaruyamaD<<<nblocks, nthreads, 0, stream>>>(pos,
										     dir,
										     indexIter,
										     d_linearV,
										     d_angularV,
										     numberParticles,
										     dt);
      CudaCheckError();
    }
    
    void forwardTime() override {
      
      if(firstStep){
	//Load all interactors into magnetic
	for(auto& interactor : this->getInteractors()){
	  magnetic->addInteractor(interactor);
	}
	
	//Load all updatables into magnetic
	for(auto& updatable : this->getUpdatables()){
	  magnetic->addUpdatable(updatable);
	}
	
	firstStep = false;
      }

      //There might be interactions that only cause force on the particles but not magnetic field,
      //and they will be ignored if we compute the forces with magneticfield = true because they
      //will not have a ForceTorqueMagneticField transverser.
      bool computeMagneticField = false;
      this->updateForce(computeMagneticField);
      computeMagneticField = true;
      this->updateForce(computeMagneticField);
      magnetic->updateMagnetization();
      updatePosition();
      this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
      magnetic->resetMagneticField();
      this->resetForce();
      this->resetTorque();
    }

  };
  
}}}}

#endif
