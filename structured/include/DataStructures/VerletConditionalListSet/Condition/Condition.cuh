#pragma once

#include "DataStructures/ExclusionsList/ExclusionsList.cuh"

namespace uammd{
namespace structured{
namespace conditions{

    class excludedConditionBase{

        protected:

            std::shared_ptr<GlobalData>            gd;
            std::shared_ptr<ExtendedParticleData>  pd;

            ////////////////////////////////////////////

            std::shared_ptr<Exclusions> exclusionList;

        public:

            ////////////////////////////////////////////

            excludedConditionBase(std::shared_ptr<GlobalData>            gd,
                                  std::shared_ptr<ExtendedParticleData>  pd,
                                  DataEntry& dataEntry):gd(gd),pd(pd){

                System::log<System::MESSAGE>("[Condition] Condition \"excludedConditionBase\" initialized");

                exclusionList = std::make_shared<Exclusions>(gd,pd,dataEntry);

            }

            ////////////////////////////////////////////

            size_t getSharedSize(){
                return exclusionList->getSharedSize();
            }
    };

}}}
