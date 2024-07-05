#include "Interactor/InteractorLoader.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    bool isInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        return Interactor::InteractorFactory::getInstance().isInteractorRegistered(potType,potSubType);
    }


    std::shared_ptr<typename uammd::Interactor>
    loadInteractor(std::shared_ptr<ExtendedSystem> sys,
                   std::shared_ptr<GlobalData>     gd,
                   std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                   std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                   std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        if(Interactor::InteractorFactory::getInstance().isInteractorNonBonded(potType,potSubType)){

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            return Interactor::InteractorFactory::getInstance().createInteractor<
                   Interactor::InteractorFactory::NonBondedCreator>(potType,potSubType,
                                                                    gd,pg,
                                                                    nl,
                                                                    data,
                                                                    path.back());
        }

        if(Interactor::InteractorFactory::getInstance().isInteractorPatchyParticle(potType,potSubType)){
            return Interactor::InteractorFactory::getInstance().createInteractor<
                   Interactor::InteractorFactory::PatchyParticleCreator
                   >(potType,potSubType,
                     gd,pg,
                     path,
                     path.back());
        }


        return Interactor::InteractorFactory::getInstance().createInteractor<
               Interactor::InteractorFactory::Creator
               >(potType,potSubType,
                 gd,pg,
                 data,
                 path.back());

    }

}}}

