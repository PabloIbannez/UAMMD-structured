#include "Utils/String/StringUtils.cuh"

namespace uammd{
namespace structured{
namespace StringUtils{

    std::string strip(const std::string& str_in){
        std::string str = str_in;

        if (str.length() == 0) {
            return str;
        }

        auto start_it = str.begin();
        auto end_it = str.rbegin();
        while (std::isspace(*start_it)) {
            ++start_it;
            if (start_it == str.end()) break;
        }
        while (std::isspace(*end_it)) {
            ++end_it;
            if (end_it == str.rend()) break;
        }
        int start_pos = start_it - str.begin();
        int end_pos = end_it.base() - str.begin();
        str = start_pos <= end_pos ? std::string(start_it, end_it.base()) : "";

        return str;
    }

}}}
