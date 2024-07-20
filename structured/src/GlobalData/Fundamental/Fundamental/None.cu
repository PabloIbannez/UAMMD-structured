#include "GlobalData/Fundamental/FundamentalHandler.cuh"
#include "GlobalData/Fundamental/FundamentalFactory.cuh"

namespace uammd{
namespace structured{
namespace Fundamental{

    class None: public FundamentalHandler{

        public:

            None(DataEntry& data):FundamentalHandler(data){}

            void updateDataEntry(DataEntry data){}

    };

}}}

REGISTER_FUNDAMENTAL(
    Fundamental,None,
    uammd::structured::Fundamental::None
)
