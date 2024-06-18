#include "GlobalData/Fundamental/FundamentalHandler.cuh"
#include "GlobalData/Fundamental/Fundamental/None.cuh"

namespace uammd{
namespace structured{
namespace Fundamental{

    None::None(DataEntry &data) : FundamentalHandler(data) {}

    void None::updateDataEntry(DataEntry data) {}

}}}
