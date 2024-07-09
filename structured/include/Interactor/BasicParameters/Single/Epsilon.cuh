#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicParameters{

    namespace Single{

        struct Epsilon{

            struct InputSingleParameters{

                std::string name;

                real epsilon;
            };

            struct SingleParameters{
                real epsilon;
            };

            template<typename T>
            static inline __host__ InputSingleParameters readSingleParameters(std::map<std::string,T>& info){

                InputSingleParameters param;

                param.name    = std::string(info.at("name"));
                param.epsilon = real(info.at("epsilon"));

                return param;
            }

            static inline __host__ SingleParameters processSingleParameters(InputSingleParameters in_par){
                SingleParameters param;

                param.epsilon = in_par.epsilon;

                return param;

            }
        };

    }

}}}}
