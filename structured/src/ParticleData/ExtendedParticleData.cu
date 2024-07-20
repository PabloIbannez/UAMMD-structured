#include "ParticleData/ExtendedParticleData.cuh"

namespace uammd{
namespace structured{

    ExtendedParticleData::ExtendedParticleData(std::shared_ptr<ExtendedSystem> sys,
                                               std::vector<std::string>       path):
        uammd::ParticleData((sys->getInput()->getDataEntry(path)).getDataSize(),sys),path(path){

        auto data = this->getSystem()->getInput()->getDataEntry(path);
        stateLoader(this, data);

    }

    ExtendedParticleData::ExtendedParticleData(std::shared_ptr<ExtendedSystem> sys):
        ExtendedParticleData(sys,{"state"}){
    }

    std::shared_ptr<ExtendedSystem> ExtendedParticleData::getSystem(){
        return getExtendedSystem(uammd::ParticleData::getSystem());
    }

    void ExtendedParticleData::updateInputState(){
        auto data = this->getSystem()->getInput()->getDataEntry(path);
        updateState(this, data);
    }

    // ExtendedParticleData cast

    std::shared_ptr<ExtendedParticleData> getExtendedParticleData(std::shared_ptr<uammd::ParticleData> pd){
        return std::static_pointer_cast<ExtendedParticleData>(pd);
    }

    std::shared_ptr<ExtendedParticleData> getExtendedParticleData(std::shared_ptr<uammd::ParticleGroup> pg){
        return std::static_pointer_cast<ExtendedParticleData>(pg->getParticleData());
    }


}}

