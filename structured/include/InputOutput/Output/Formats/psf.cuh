#ifndef __PSF__
#define __PSF__

namespace uammd{
namespace structured{
namespace psf{

    void WritePSF(std::shared_ptr<ParticleGroup> pg,
                  std::ofstream& out){

        auto pd  = pg->getParticleData();
        auto sys = pd->getSystem();

        int N = pg->getNumberParticles();

        auto id    = pd->getId(access::location::cpu,
                               access::mode::read);

        auto radius  = pd->getRadius(access::location::cpu,
                                     access::mode::read);
        auto mass    = pd->getMass(access::location::cpu,
                                   access::mode::read);
        auto charge  = pd->getCharge(access::location::cpu,
                                     access::mode::read);

        auto pos   = pd->getPos(access::location::cpu,
                                access::mode::read);

        auto resId   = pd->getResId(access::location::cpu,
                                    access::mode::read);
        auto chainId = pd->getChainId(access::location::cpu,
                                      access::mode::read);
        auto modelId = pd->getModelId(access::location::cpu,
                                      access::mode::read);

        auto batId = pd->getBatchId(access::location::cpu,
                                    access::mode::read);

        auto groupIndex  = pg->getIndexIterator(access::location::cpu);
        auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

        std::map<std::tuple<int,int,int>,int> segid;

        int segidCount = 0;
        fori(0,N){
            std::tuple<int,int,int> mdlCh = std::make_tuple(batId[i],modelId[i],chainId[i]);

            if(segid.count(mdlCh) == 0){
                segid[mdlCh] = segidCount;
                segidCount++;
            }
        }

        for(auto si : segid){
            System::log<System::DEBUG>("[OutPut PSF] (Bat,Model,Chain) SegId: (%i,%i,%i) %i",
                                      std::get<0>(si.first),
                                      std::get<1>(si.first),
                                      std::get<2>(si.first),si.second);
        }

        std::map<int,real> type2radius;
        fori(0,N){
            type2radius[int(pos[i].w)] = radius[i];
        }

        out << "PSF CMAP" << std::endl << std::endl;

        out << "       " << type2radius.size() << " !NTITLE" << std::endl;

        for(auto& tr : type2radius){
            out << "REMARKS " << tr.first << " " << tr.second << std::endl;
        }

        out << std::endl;

        out << std::fixed              <<
               std::right              <<
               std::setw(8) << N       << " !NATOM" << std::endl;

        std::map<int,int> id_index;
        fori(0,N){
            //int id_   = id[groupIndex[i]];
            int id_   = i;
            int index = sortedIndex[id[groupIndex[i]]];

            id_index[id_+1]=index; //Starts at 1 !!! psf format stuff
        }

        for(const auto& ii : id_index){

            int index = ii.second;
            int id_   = ii.first;

            std::tuple<int,int,int> batMdlCh = std::make_tuple(batId[index],modelId[index],chainId[index]);

            out << std::fixed            <<
                   std::right            <<
                   std::setw(8) << id_   <<
                   " "                   <<
                   std::fixed            <<
                   std::left             <<
                   std::setw(4) << segid[batMdlCh] <<
                   " "                   <<
                   std::setw(4) << modelId[index] <<
                   " "                   <<
                   std::setw(4) << resId[index] <<
                   " "                   <<
                   //std::setw(4) << typeParameterHandler->getTypeParameters(int(pos[index].w)).name <<
                   std::setw(4) << int(pos[index].w) <<
                   " "                   <<
                   //std::setw(4) << typeParameterHandler->getTypeParameters(int(pos[index].w)).name <<
                   std::setw(4) << int(pos[index].w) <<
                   " "                   <<
                   //std::setw(4) << int(pos[index].w) <<
                   //" "                   <<
                   std::fixed            <<
                   std::right            <<
                   std::setprecision(6)  <<
                   std::setw(14) << charge[index] <<
                   std::setprecision(6)  <<
                   std::setw(14) << mass[index] <<
                   "       0"  << std::endl;
        }

        std::vector<int> ids;
        for(const auto& ii : id_index){
            ids.push_back(ii.first);
        }

        std::vector<std::pair<int,int>> bonds;
        for(uint i=0;i<ids.size()-1;i++){
            int id_     = ids[i];
            int id_next = ids[i+1];

            int index      = id_index[id_];
            int index_next = id_index[id_next];

            std::tuple<int,int,int> batMdlCh      = std::make_tuple(batId[index],modelId[index],chainId[index]);
            std::tuple<int,int,int> batMdlCh_next = std::make_tuple(batId[index_next],modelId[index_next],chainId[index_next]);

            if(batMdlCh == batMdlCh_next and (resId[index]+1)==resId[index_next]){
                bonds.push_back(std::make_pair(id_,id_next));
            }
        }

        out << std::endl;
        out << std::endl;

        out << std::fixed                   <<
               std::right                   <<
               std::setw(10) << bonds.size() << " !NBOND: bonds" << std::endl;

        int i=0;
        for(const auto& bond : bonds){

            i++;
            out << std::fixed                   <<
                   std::right                   <<
                   std::setw(10) << bond.first  <<
                   std::fixed                   <<
                   std::right                   <<
                   std::setw(10) << bond.second ;

            if(i%4==0){
                out << std::endl;
            }
        }

        out << std::endl;
        out << std::endl;
    }

}}}

#endif
