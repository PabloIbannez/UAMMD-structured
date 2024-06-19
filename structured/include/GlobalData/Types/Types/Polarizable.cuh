#ifndef __POLARIZABLE_TYPE__
#define __POLARIZABLE_TYPE__

namespace uammd{
namespace structured{
namespace Types{

    struct Polarizable_{
        template<typename T>
        static void loadType(std::map<std::string,std::map<std::string,real>>& nameToData,
                             std::map<std::string,T>& typeData);

        static void loadTypesIntoParticleData(std::shared_ptr<ParticleData> pd,
                                              std::map<int,std::string>&    idToName,
                                              std::map<std::string,std::map<std::string,real>>& nameToData);
    };

    using Polarizable = Types_<Polarizable_>;

}}}

#endif
