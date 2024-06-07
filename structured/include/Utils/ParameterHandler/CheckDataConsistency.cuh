#ifndef CHECK_DATA_CONSISTENCY_CUH
#define CHECK_DATA_CONSISTENCY_CUH

namespace uammd{
namespace structured{
namespace Batching{

    bool checkDataConsistency(DataEntry& data){

        //Check "batchId" is consistent
        std::vector<int> batchIds = data.getData<int>("batchId");

        //Convert to set and check is iota
        std::set<int>    batchIdsSet(batchIds.begin(), batchIds.end());
        std::vector<int> batchIdsSetVector(batchIdsSet.begin(), batchIdsSet.end());

        std::vector<int> batchIdsIota(batchIdsSet.size());
        std::iota(batchIdsIota.begin(), batchIdsIota.end(), 0);

        if(batchIdsSetVector != batchIdsIota){
            std::vector<int> diff(batchIdsSetVector.size() + batchIdsIota.size());
            auto it = std::set_symmetric_difference(batchIdsSetVector.begin(), batchIdsSetVector.end(),
                                                    batchIdsIota.begin(), batchIdsIota.end(), diff.begin());
            diff.resize(it - diff.begin());

            std::string msg = "";
            for(auto i : diff){
                msg += std::to_string(i) + " ";
            }
            System::log<System::WARNING>("[Batching] batchId is not iota. Missing: %s", msg.c_str());
            return false;
        }

        //Check there are the same number of different "batchId"
        //Count the number of entries with each "batchId" is equal to zero
        int nBatchId = std::count(batchIds.begin(), batchIds.end(), 0);
        for(auto it = batchIds.begin(); it != batchIds.end(); ++it){
            if(std::count(batchIds.begin(), batchIds.end(), *it) != nBatchId){
                System::log<System::WARNING>("[Batching] There are different number of entries with the same batchId");
                return false;
            }
        }


        return true;
    }

    bool checkParticleGroupDataConsistency(std::shared_ptr<ParticleGroup> pg,
                                           DataEntry& data){

        auto pd = pg->getParticleData();
        auto pgBatchId = pg->getPropertyIterator(pd->getBatchId(access::location::cpu, access::mode::read).begin(),
                                                 access::location::cpu);

        std::vector<int> dataBatchId = data.getData<int>("batchId");

        //

        std::set<int> batchIdsSet;
        for(int i = 0; i < pg->getNumberParticles(); ++i){
            batchIdsSet.insert(pgBatchId[i]);
        }
        std::set<int> dataBatchIdSet(dataBatchId.begin(), dataBatchId.end());

        //Check if the two sets are equal (order does not matter)
        if(batchIdsSet != dataBatchIdSet){
            System::log<System::WARNING>("[Batching] batchId is not the same in ParticleGroup and DataEntry");
            return false;
        }

        return true;
    }

}}}

#endif
