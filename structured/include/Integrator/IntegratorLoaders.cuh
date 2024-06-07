#ifndef __INTEGRATOR_LOADER__
#define __INTEGRATOR_LOADER__
namespace uammd{
namespace structured{
namespace IntegratorLoader{

    bool isIntegratorAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string integratorType    = data.getType();
        std::string integratorSubType = data.getSubType();
        if("SteepestDescent" == integratorType and "SteepestDescent" == integratorSubType){
            return true;
        }
        if("Brownian" == integratorType and "EulerMaruyama" == integratorSubType){
            return true;
        }
        if("Brownian" == integratorType and "EulerMaruyamaRigidBody" == integratorSubType){
            return true;
        }
        if("Brownian" == integratorType and "EulerMaruyamaRigidBodyPatchesState" == integratorSubType){
            return true;
        }
        if("Langevin" == integratorType and "BBK" == integratorSubType){
            return true;
        }
        if("Langevin" == integratorType and "GJF" == integratorSubType){
            return true;
        }
        if("DPD" == integratorType and "Verlet" == integratorSubType){
            return true;
        }
        if("SPH" == integratorType and "Verlet" == integratorSubType){
            return true;
        }
        if("BDHIOpenBoundary" == integratorType and "Cholesky" == integratorSubType){
            return true;
        }
        if("BDHIOpenBoundary" == integratorType and "Lanczos" == integratorSubType){
            return true;
        }
        if("BDHITriplyPeriodic" == integratorType and "ForceCouplingMethod" == integratorSubType){
            return true;
        }
        if("BDHITriplyPeriodic" == integratorType and "PositivelySplitEwald" == integratorSubType){
            return true;
        }
        if("BDHITriplyPeriodic" == integratorType and "FluctuatingImmersedBoundary" == integratorSubType){
            return true;
        }
        if("BDHIDoublyPeriodic" == integratorType and "DPStokes" == integratorSubType){
            return true;
        }
        if("FluctuatingHydrodynamics" == integratorType and "CompressibleInertialCoupling" == integratorSubType){
            return true;
        }
        if("FluctuatingHydrodynamics" == integratorType and "IncompressibleInertialCoupling" == integratorSubType){
            return true;
        }
        if("Magnetic" == integratorType and "Brownian" == integratorSubType){
            return true;
        }
        if("Magnetic" == integratorType and "Fixed" == integratorSubType){
            return true;
        }
        if("Magnetic" == integratorType and "ForceCouplingMethod" == integratorSubType){
            return true;
        }
        if("BDHIDoublyPeriodic" == integratorType and "Quasi2D" == integratorSubType){
            return true;
        }
        if("Verlet" == integratorType and "VelocityVerlet" == integratorSubType){
            return true;
        }
        if("None" == integratorType and "None" == integratorSubType){
            return true;
        }
        
        return false;

    }

    
    std::shared_ptr<typename uammd::Integrator>
    loadIntegrator(std::shared_ptr<ExtendedSystem> sys,
                   std::shared_ptr<GlobalData>     gd,
                   std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                   std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string integratorType    = data.getType();
        std::string integratorSubType = data.getSubType();

        std::shared_ptr<typename uammd::Integrator> integrator;
        bool found = false;
        
        if("SteepestDescent" == integratorType and "SteepestDescent" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected Minimization::SteepestDescent::SteepestDescent integrator",path.back().c_str());
            integrator = std::make_shared<Minimization::SteepestDescent::SteepestDescent>(gd,pg,data,path.back());
            found = true;
        }
        if("Brownian" == integratorType and "EulerMaruyama" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Brownian::EulerMaruyama integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Brownian::EulerMaruyama>(gd,pg,data,path.back());
            found = true;
        }
        if("Brownian" == integratorType and "EulerMaruyamaRigidBody" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Brownian::EulerMaruyamaRigidBody integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Brownian::EulerMaruyamaRigidBody>(gd,pg,data,path.back());
            found = true;
        }
        if("Brownian" == integratorType and "EulerMaruyamaRigidBodyPatchesState" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Brownian::EulerMaruyamaRigidBodyPatchesState integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Brownian::EulerMaruyamaRigidBodyPatchesState>(gd,pg,data,path.back());
            found = true;
        }
        if("Langevin" == integratorType and "BBK" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Langevin::BBK integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Langevin::BBK>(gd,pg,data,path.back());
            found = true;
        }
        if("Langevin" == integratorType and "GJF" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Langevin::GJF integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Langevin::GJF>(gd,pg,data,path.back());
            found = true;
        }
        if("DPD" == integratorType and "Verlet" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::DPD::Verlet integrator",path.back().c_str());
            integrator = std::make_shared<NVT::DPD::Verlet>(gd,pg,data,path.back());
            found = true;
        }
        if("SPH" == integratorType and "Verlet" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::SPH::Verlet integrator",path.back().c_str());
            integrator = std::make_shared<NVT::SPH::Verlet>(gd,pg,data,path.back());
            found = true;
        }
        if("BDHIOpenBoundary" == integratorType and "Cholesky" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::BDHIOpenBoundary::Cholesky integrator",path.back().c_str());
            integrator = std::make_shared<NVT::BDHIOpenBoundary::Cholesky>(gd,pg,data,path.back());
            found = true;
        }
        if("BDHIOpenBoundary" == integratorType and "Lanczos" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::BDHIOpenBoundary::Lanczos integrator",path.back().c_str());
            integrator = std::make_shared<NVT::BDHIOpenBoundary::Lanczos>(gd,pg,data,path.back());
            found = true;
        }
        if("BDHITriplyPeriodic" == integratorType and "ForceCouplingMethod" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::BDHITriplyPeriodic::ForceCouplingMethod integrator",path.back().c_str());
            integrator = std::make_shared<NVT::BDHITriplyPeriodic::ForceCouplingMethod>(gd,pg,data,path.back());
            found = true;
        }
        if("BDHITriplyPeriodic" == integratorType and "PositivelySplitEwald" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::BDHITriplyPeriodic::PositivelySplitEwald integrator",path.back().c_str());
            integrator = std::make_shared<NVT::BDHITriplyPeriodic::PositivelySplitEwald>(gd,pg,data,path.back());
            found = true;
        }
        if("BDHITriplyPeriodic" == integratorType and "FluctuatingImmersedBoundary" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::BDHITriplyPeriodic::FluctuatingImmersedBoundary integrator",path.back().c_str());
            integrator = std::make_shared<NVT::BDHITriplyPeriodic::FluctuatingImmersedBoundary>(gd,pg,data,path.back());
            found = true;
        }
        if("BDHIDoublyPeriodic" == integratorType and "DPStokes" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::BDHIDoublyPeriodic::DPStokes integrator",path.back().c_str());
            integrator = std::make_shared<NVT::BDHIDoublyPeriodic::DPStokes>(gd,pg,data,path.back());
            found = true;
        }
        if("FluctuatingHydrodynamics" == integratorType and "CompressibleInertialCoupling" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::FluctuatingHydrodynamics::CompressibleInertialCoupling integrator",path.back().c_str());
            integrator = std::make_shared<NVT::FluctuatingHydrodynamics::CompressibleInertialCoupling>(gd,pg,data,path.back());
            found = true;
        }
        if("FluctuatingHydrodynamics" == integratorType and "IncompressibleInertialCoupling" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::FluctuatingHydrodynamics::IncompressibleInertialCoupling integrator",path.back().c_str());
            integrator = std::make_shared<NVT::FluctuatingHydrodynamics::IncompressibleInertialCoupling>(gd,pg,data,path.back());
            found = true;
        }
        if("Magnetic" == integratorType and "Brownian" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Magnetic::Brownian integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Magnetic::Brownian>(gd,pg,data,path.back());
            found = true;
        }
        if("Magnetic" == integratorType and "Fixed" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Magnetic::Fixed integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Magnetic::Fixed>(gd,pg,data,path.back());
            found = true;
        }
        if("Magnetic" == integratorType and "ForceCouplingMethod" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::Magnetic::ForceCouplingMethod integrator",path.back().c_str());
            integrator = std::make_shared<NVT::Magnetic::ForceCouplingMethod>(gd,pg,data,path.back());
            found = true;
        }
        if("BDHIDoublyPeriodic" == integratorType and "Quasi2D" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVT::BDHIDoublyPeriodic::Quasi2D integrator",path.back().c_str());
            integrator = std::make_shared<NVT::BDHIDoublyPeriodic::Quasi2D>(gd,pg,data,path.back());
            found = true;
        }
        if("Verlet" == integratorType and "VelocityVerlet" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected NVE::Verlet::VelocityVerlet integrator",path.back().c_str());
            integrator = std::make_shared<NVE::Verlet::VelocityVerlet>(gd,pg,data,path.back());
            found = true;
        }
        if("None" == integratorType and "None" == integratorSubType){
            System::log<System::MESSAGE>("[IntegratorLoader] (%s) Detected Special::None::None integrator",path.back().c_str());
            integrator = std::make_shared<Special::None::None>(gd,pg,data,path.back());
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[IntegratorLoader] (%s) Could not find integrator %s::%s",
                                            path.back().c_str(),integratorType.c_str(),integratorSubType.c_str());
        }

        return integrator;

    }

    }}}
#endif
