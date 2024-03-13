#ifndef __EXTENDED_PARTICLE_DATA__
#define __EXTENDED_PARTICLE_DATA__

namespace uammd{
namespace structured{

    class ExtendedParticleData: public uammd::ParticleData{

        private:

            std::vector<std::string> path;

        public:

            ExtendedParticleData(std::shared_ptr<ExtendedSystem> sys,
                                 std::vector<std::string>       path):uammd::ParticleData((sys->getInput()->getDataEntry(path)).getDataSize(),sys),
                                                                      path(path){

                auto data = this->getSystem()->getInput()->getDataEntry(path);
                stateLoader(this, data);

            }

            ExtendedParticleData(std::shared_ptr<ExtendedSystem> sys):ExtendedParticleData(sys,{"state"}){
            }

            std::shared_ptr<ExtendedSystem> getSystem(){
                return getExtendedSystem(uammd::ParticleData::getSystem());
            }

            void updateInputState(){
                auto data = this->getSystem()->getInput()->getDataEntry(path);
                updateState(this, data);
            }
    };

    std::shared_ptr<ExtendedParticleData> getExtendedParticleData(std::shared_ptr<uammd::ParticleData> pd){
        return std::static_pointer_cast<ExtendedParticleData>(pd);
    }

    std::shared_ptr<ExtendedParticleData> getExtendedParticleData(std::shared_ptr<uammd::ParticleGroup> pg){
        return std::static_pointer_cast<ExtendedParticleData>(pg->getParticleData());
    }


}}

#endif

