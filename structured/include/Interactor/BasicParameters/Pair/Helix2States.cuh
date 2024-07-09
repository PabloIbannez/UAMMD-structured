#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace BasicParameters{

    namespace Pairs{

        template<typename potential>
        struct Helix2States {

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real Eb0;
                typename potential::params helixParams0;

                tensor3 R_H0;

                real Eb1;
                typename potential::params helixParams1;

                tensor3 R_H1;

                real prob_0_to_1;
                real prob_1_to_0;
            };

            struct PairParameters{

                real Eb0;
                typename potential::params helixParams0;

                tensor3 R_H0;

                real Eb1;
                typename potential::params helixParams1;

                tensor3 R_H1;

                real prob_0_to_1;
                real prob_1_to_0;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                //State 0

                param.Eb0  = real(info.at("Eb0"));

                param.helixParams0 = potential::readParamsMap(info,"0");

                real3 e_x0 = real3(info.at("e_x0"));
                real3 e_y0 = real3(info.at("e_y0"));
                real3 e_z0 = real3(info.at("e_z0"));

                param.R_H0 = Helix<potential>::loadR_H(e_x0,e_y0,e_z0);

                //State 1

                param.Eb1  = real(info.at("Eb1"));

                param.helixParams1 = potential::readParamsMap(info,"1");

                real3 e_x1 = real3(info.at("e_x1"));
                real3 e_y1 = real3(info.at("e_y1"));
                real3 e_z1 = real3(info.at("e_z1"));

                param.R_H1 = Helix<potential>::loadR_H(e_x1,e_y1,e_z1);

                //Probabilities

                param.prob_0_to_1 = real(info.at("prob_0_to_1"));
                param.prob_1_to_0 = real(info.at("prob_1_to_0"));

                //Print

                System::log<System::MESSAGE>("[HelixParameters] Added parameter Eb0: %f",param.Eb0);

                System::log<System::MESSAGE>("               %f,%f,%f",param.R_H0.xx,param.R_H0.xy,param.R_H0.xz);
                System::log<System::MESSAGE>("[Helix] R_H0 = %f,%f,%f",param.R_H0.yx,param.R_H0.yy,param.R_H0.yz);
                System::log<System::MESSAGE>("               %f,%f,%f",param.R_H0.zx,param.R_H0.zy,param.R_H0.zz);

                ////////////////////////////////////////////////////////////////////////////////////

                System::log<System::MESSAGE>("[HelixParameters] Added parameter Eb1: %f",param.Eb1);

                System::log<System::MESSAGE>("               %f,%f,%f",param.R_H1.xx,param.R_H1.xy,param.R_H1.xz);
                System::log<System::MESSAGE>("[Helix] R_H1 = %f,%f,%f",param.R_H1.yx,param.R_H1.yy,param.R_H1.yz);
                System::log<System::MESSAGE>("               %f,%f,%f",param.R_H1.zx,param.R_H1.zy,param.R_H1.zz);

                ////////////////////////////////////////////////////////////////////////////////////

                System::log<System::MESSAGE>("[HelixParameters] Added parameter prob_0_to_1: %f",param.prob_0_to_1);
                System::log<System::MESSAGE>("[HelixParameters] Added parameter prob_1_to_0: %f",param.prob_1_to_0);

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                //State 0

                param.Eb0  = in_par.Eb0;
                param.helixParams0 = in_par.helixParams0;
                param.R_H0 = in_par.R_H0;

                //State 1

                param.Eb1  = in_par.Eb1;
                param.helixParams1 = in_par.helixParams1;
                param.R_H1 = in_par.R_H1;

                //Probabilities

                param.prob_0_to_1 = in_par.prob_0_to_1;
                param.prob_1_to_0 = in_par.prob_1_to_0;

                return param;
            }
        };
    }

}}}}
