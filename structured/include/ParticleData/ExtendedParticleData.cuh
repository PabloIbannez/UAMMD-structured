#pragma once

#include"uammd.cuh"

#include"System/ExtendedSystem.cuh"

#include"ParticleData/ParticleData.cuh"
#include"ParticleData/ParticleGroup.cuh"
#include"ParticleData/StateLoader.cuh"

#include<vector>
#include<string>
#include<memory>

namespace uammd{
namespace structured{

    class ExtendedParticleData: public uammd::ParticleData{

        private:

            std::vector<std::string> path;

        public:

            ExtendedParticleData(std::shared_ptr<ExtendedSystem> sys,
                                 std::vector<std::string>       path);

            ExtendedParticleData(std::shared_ptr<ExtendedSystem> sys);

            std::shared_ptr<ExtendedSystem> getSystem();

            void updateInputState();
    };

    std::shared_ptr<ExtendedParticleData> getExtendedParticleData(std::shared_ptr<uammd::ParticleData> pd);
    std::shared_ptr<ExtendedParticleData> getExtendedParticleData(std::shared_ptr<uammd::ParticleGroup> pg);

}}

