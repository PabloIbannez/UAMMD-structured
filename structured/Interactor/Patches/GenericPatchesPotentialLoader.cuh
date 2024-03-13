#ifndef __GENERIC_PATCHY_PARTICLES_POTENTIALS_LOADER__
#define __GENERIC_PATCHY_PARTICLES_POTENTIALS_LOADER__
namespace uammd{
namespace structured{
namespace Potentials{
namespace GenericPatchyParticlesLoader{

    bool isPatchyParticlesInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                                              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();
        
        if("SurfacePatches" == potType and "Linker" == potSubType){
            return true;
        }
            
        if("NonBondedPatches" == potType and "HelixExponential" == potSubType){
            return true;
        }
            
        if("NonBondedPatches" == potType and "HelixExponential2States" == potSubType){
            return true;
        }
            
        if("NonBondedPatches" == potType and "HelixCosine" == potSubType){
            return true;
        }
            
        if("NonBondedPatches" == potType and "HelixCosine2States" == potSubType){
            return true;
        }
            
        if("NonBondedPatches" == potType and "DistanceSwitchExponential" == potSubType){
            return true;
        }
            
        if("NonBondedPatches" == potType and "DistanceSwitchCosine" == potSubType){
            return true;
        }
            
        return false;

    }

    
    std::shared_ptr<typename uammd::Interactor>
    loadGenericPatchyParticles(std::shared_ptr<ExtendedSystem> sys,
                               std::shared_ptr<GlobalData>  gd,std::shared_ptr<ParticleGroup>  pg,
                               std::shared_ptr<GlobalData>  patchesGd,std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                               std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> patchesPg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        std::shared_ptr<typename uammd::Interactor> interactor;
        bool found = false;

        if("NonBondedPatches" == potType and "HelixExponential" == potSubType){
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected NonBondedPatches::HelixExponential potential",path.back().c_str());
            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);
            std::shared_ptr<Potentials::NonBondedPatches::HelixExponential> pot = std::make_shared<Potentials::NonBondedPatches::HelixExponential>(gd,pg,patchesGd,patchesPg,data);
            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));
            interactor = std::make_shared<typename Interactor::PairInteractor<Potentials::NonBondedPatches::HelixExponential,VerletConditionalListSetBase>>(patchesGd,patchesPg,data,pot,nl,path.back());
            found = true;
        }
        if("NonBondedPatches" == potType and "HelixExponential2States" == potSubType){
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected NonBondedPatches::HelixExponential2States potential",path.back().c_str());
            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);
            std::shared_ptr<Potentials::NonBondedPatches::HelixExponential2States> pot = std::make_shared<Potentials::NonBondedPatches::HelixExponential2States>(gd,pg,patchesGd,patchesPg,data);
            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));
            interactor = std::make_shared<typename Interactor::PairInteractor<Potentials::NonBondedPatches::HelixExponential2States,VerletConditionalListSetBase>>(patchesGd,patchesPg,data,pot,nl,path.back());
            found = true;
        }
        if("NonBondedPatches" == potType and "HelixCosine" == potSubType){
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected NonBondedPatches::HelixCosine potential",path.back().c_str());
            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);
            std::shared_ptr<Potentials::NonBondedPatches::HelixCosine> pot = std::make_shared<Potentials::NonBondedPatches::HelixCosine>(gd,pg,patchesGd,patchesPg,data);
            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));
            interactor = std::make_shared<typename Interactor::PairInteractor<Potentials::NonBondedPatches::HelixCosine,VerletConditionalListSetBase>>(patchesGd,patchesPg,data,pot,nl,path.back());
            found = true;
        }
        if("NonBondedPatches" == potType and "HelixCosine2States" == potSubType){
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected NonBondedPatches::HelixCosine2States potential",path.back().c_str());
            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);
            std::shared_ptr<Potentials::NonBondedPatches::HelixCosine2States> pot = std::make_shared<Potentials::NonBondedPatches::HelixCosine2States>(gd,pg,patchesGd,patchesPg,data);
            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));
            interactor = std::make_shared<typename Interactor::PairInteractor<Potentials::NonBondedPatches::HelixCosine2States,VerletConditionalListSetBase>>(patchesGd,patchesPg,data,pot,nl,path.back());
            found = true;
        }
        if("NonBondedPatches" == potType and "DistanceSwitchExponential" == potSubType){
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected NonBondedPatches::DistanceSwitchExponential potential",path.back().c_str());
            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);
            std::shared_ptr<Potentials::NonBondedPatches::DistanceSwitchExponential> pot = std::make_shared<Potentials::NonBondedPatches::DistanceSwitchExponential>(gd,pg,patchesGd,patchesPg,data);
            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));
            interactor = std::make_shared<typename Interactor::PairInteractor<Potentials::NonBondedPatches::DistanceSwitchExponential,VerletConditionalListSetBase>>(patchesGd,patchesPg,data,pot,nl,path.back());
            found = true;
        }
        if("NonBondedPatches" == potType and "DistanceSwitchCosine" == potSubType){
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected NonBondedPatches::DistanceSwitchCosine potential",path.back().c_str());
            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);
            std::shared_ptr<Potentials::NonBondedPatches::DistanceSwitchCosine> pot = std::make_shared<Potentials::NonBondedPatches::DistanceSwitchCosine>(gd,pg,patchesGd,patchesPg,data);
            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));
            interactor = std::make_shared<typename Interactor::PairInteractor<Potentials::NonBondedPatches::DistanceSwitchCosine,VerletConditionalListSetBase>>(patchesGd,patchesPg,data,pot,nl,path.back());
            found = true;
        }

        if("SurfacePatches" == potType and "Linker" == potSubType){
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected SurfacePatches::Linker potential",path.back().c_str());
            std::shared_ptr<Potentials::SurfacePatches::Linker> pot = std::make_shared<Potentials::SurfacePatches::Linker>(gd,pg,patchesGd,patchesPg,data);
            interactor = std::make_shared<typename Interactor::SingleInteractor<Potentials::SurfacePatches::Linker>>(patchesGd,patchesPg,data,pot,path.back());
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[GenericPatchyParticles] (%s) Potential of type %s and subType %s not found",
                                           path.back().c_str(),data.getType().c_str(),data.getSubType().c_str());
        }

        return interactor;

    }

    }}}}
#endif
