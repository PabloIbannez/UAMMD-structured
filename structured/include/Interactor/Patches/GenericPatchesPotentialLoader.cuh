#ifndef __GENERIC_PATCHY_PARTICLES_POTENTIALS_LOADER__
#define __GENERIC_PATCHY_PARTICLES_POTENTIALS_LOADER__
namespace uammd{
namespace structured{
namespace Potentials{
namespace GenericPatchyParticlesLoader{

    inline
    bool isPatchyParticlesInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                                              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();


        return false;

    }


    inline
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

        if(not found){
            System::log<System::CRITICAL>("[GenericPatchyParticles] (%s) Potential of type %s and subType %s not found",
                                           path.back().c_str(),data.getType().c_str(),data.getSubType().c_str());
        }

        return interactor;

    }

    }}}}
#endif
