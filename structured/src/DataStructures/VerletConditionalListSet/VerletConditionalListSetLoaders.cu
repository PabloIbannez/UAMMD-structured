#include "DataStructures/VerletConditionalListSet/VerletConditionalListSetLoaders.cuh"

namespace uammd{
namespace structured{
namespace VerletConditionalListSetLoaders{

    std::shared_ptr<uammd::structured::VerletConditionalListSetBase>
    loadVerletConditionalListSet(std::shared_ptr<ExtendedSystem> sys,
                                 std::shared_ptr<GlobalData>     gd,
                                 std::map<std::string,std::shared_ptr<ParticleGroup>> groups,
                                 std::vector<std::string>       path){


        DataEntry data = sys->getInput()->getDataEntry(path);
        //Check data type is VerletConditionalListSet
        if(data.getType() != "VerletConditionalListSet"){
            System::log<System::CRITICAL>("[VerletConditionalListSetLoader] (%s) Data type is not VerletConditionalListSet!",path.back().c_str());
        }

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        //////////////////////////////////////////////////////////////////////////////

        std::string condition = data.getSubType();

        return VerletConditionalListSetFactory::getInstance().createVerletConditionalListSet(condition,
                                                                                             gd,pg,
                                                                                             data,
                                                                                             path.back());

    }
}}}
