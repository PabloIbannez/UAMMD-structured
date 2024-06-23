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


        if(not found){
            System::log<System::CRITICAL>("[IntegratorLoader] (%s) Could not find integrator %s::%s",
                                            path.back().c_str(),integratorType.c_str(),integratorSubType.c_str());
        }

        return integrator;

    }

    }}}
#endif
