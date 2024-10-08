#include "Integrator/IntegratorLoaders.cuh"

namespace uammd{
namespace structured{
namespace IntegratorLoader{

    bool isIntegratorAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string integratorType    = data.getType();
        std::string integratorSubType = data.getSubType();

        return Integrator::IntegratorFactory::getInstance().isIntegratorRegistered(integratorType,integratorSubType);
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

        //////////////////////////////////////////////////////////////

        return Integrator::IntegratorFactory::getInstance().createIntegrator(integratorType,integratorSubType,
                                                                 gd,pg,
                                                                 data,
                                                                 path.back());
    }

}}}
