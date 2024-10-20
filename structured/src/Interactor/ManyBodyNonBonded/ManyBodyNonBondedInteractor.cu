#pragma once

#include "Interactor/ManyBodyNonBonded/ManyBodyNonBondedInteractor.cuh"

namespace uammd {
namespace structured {
namespace Interactor {

    namespace ManyBodyNonBonded_ns{

    }

    template<int THREADS_PER_BLOCK>
    ManyBodyNonBonded_<THREADS_PER_BLOCK>::ManyBodyNonBonded_(std::shared_ptr<GlobalData>    gd,
                                                              std::shared_ptr<ParticleGroup> pg,
                                                              std::vector<std::string>     path,
                                                              std::string name):Interactor(pg,"ManyBodyNonBondedInteractor: \"" +name+"\""),
                                                                                           gd(gd){
        pd  = getExtendedParticleData(this->pg);

        // TODO ...
    }

    template<int THREADS_PER_BLOCK>
    void ManyBodyNonBonded_<THREADS_PER_BLOCK>::sum(uammd::Interactor::Computables comp,cudaStream_t st){

        // TODO ...
    }

}}}
