#include "Interactor/Patches/PatchInteractorLoader.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    bool isPatchInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                                    std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        return PatchesFactory::getInstance().isInteractorRegistered(potType,potSubType);
    }


    std::shared_ptr<typename uammd::Interactor>
    loadPatchInteractor(std::shared_ptr<ExtendedSystem> sys,
                        std::shared_ptr<GlobalData>  gd,std::shared_ptr<ParticleGroup>  pg,
                        std::shared_ptr<GlobalData>  patchesGd,std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                        std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                        std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> patchesPg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        if(PatchesFactory::getInstance().isInteractorNonBonded(potType,potSubType)){

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            return PatchesFactory::getInstance().createInteractor<
                   PatchesFactory::NonBondedPatchesCreator>(potType,potSubType,
                                                            gd,pg,
                                                            patchesGd,patchesPg,
                                                            nl,
                                                            data,
                                                            path.back());
        }

        return PatchesFactory::getInstance().createInteractor<
               PatchesFactory::PatchesCreator
               >(potType,potSubType,
                 gd,pg,
                 patchesGd,patchesPg,
                 data,
                 path.back());

    }


}}}

