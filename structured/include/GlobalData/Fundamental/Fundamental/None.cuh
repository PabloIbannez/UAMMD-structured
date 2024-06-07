#ifndef __NONE_FUNDAMENTAL__
#define __NONE_FUNDAMENTAL__

namespace uammd{
namespace structured{
namespace Fundamental{

    class None: public FundamentalHandler{

        public:

            None(DataEntry& data):FundamentalHandler(data){}

            void updateDataEntry(DataEntry data){}

    };

}}}

#endif
