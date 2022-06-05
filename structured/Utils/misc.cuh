#ifndef __MISC___
#define __MISC___

namespace uammd{
namespace structured{
namespace Miscellany{

std::vector<std::string> split(std::string str,std::string delimiter){
    
    std::vector<std::string> splitBuffer;

    size_t pos = 0;
    std::string s = str;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        splitBuffer.push_back(s.substr(0, pos));
        s.erase(0, pos + std::string(delimiter).length());
    }
    splitBuffer.push_back(s);

    return splitBuffer;
}

std::string str2str(std::string str,std::shared_ptr<System> sys){
    sys->template log<System::MESSAGE>("[String2string] "
                                       "Input string: %s, str: %s",str.c_str(),str.c_str());
    return str;
}

real str2real(std::string str,std::shared_ptr<System> sys){
    real f = std::stof(str);
    sys->template log<System::MESSAGE>("[String2real] "
                                       "Input string: %s, real: %f",str.c_str(),f);
    return f;
}

real3 str2real3(std::string str,std::shared_ptr<System> sys){
    real3 f3 = make_real3(0);
    sys->template log<System::MESSAGE>("[String2real3] "
                                       "Input string: %s, real3: %f %f %f",str.c_str(),f3.x,f3.y,f3.z);
    return f3;
}

int str2int(std::string str,std::shared_ptr<System> sys){
    int i = std::stoi(str);
    sys->template log<System::MESSAGE>("[String2real] "
                                       "Input string: %s, int: %d",str.c_str(),i);
    return i;
}

}}}

#endif
