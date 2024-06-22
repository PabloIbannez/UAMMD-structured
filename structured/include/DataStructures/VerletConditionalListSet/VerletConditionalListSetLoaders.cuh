#ifndef __VERLET_CONDITIONAL_LIST_SET_LOADERS__
#define __VERLET_CONDITIONAL_LIST_SET_LOADERS__

namespace uammd{
namespace structured{
namespace VerletConditionalListSetLoaders{

    std::shared_ptr<uammd::structured::VerletConditionalListSetBase>
    loadVerletConditionalListSet(std::shared_ptr<ExtendedSystem> sys,
                                 std::shared_ptr<GlobalData>    gd,
                                 std::map<std::string,std::shared_ptr<ParticleGroup>> groups,
                                 std::vector<std::string>       path);
}}}
#endif
