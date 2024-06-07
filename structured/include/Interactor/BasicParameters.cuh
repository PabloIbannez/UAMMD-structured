#ifndef __BASIC_PARAMETERS__
#define __BASIC_PARAMETERS__

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

        struct LennardJones{

            struct InputSingleParameters{

                std::string name;

                real epsilon;
                real sigma;
            };

            struct SingleParameters{
                real epsilon;
                real sigma;
            };

            template<typename T>
            static inline __host__ InputSingleParameters readSingleParameters(std::map<std::string,T>& info){

                InputSingleParameters param;

                param.name    = std::string(info.at("name"));
                param.epsilon = real(info.at("epsilon"));
                param.sigma   = real(info.at("sigma"));

                return param;
            }

            static inline __host__ SingleParameters processSingleParameters(InputSingleParameters in_par){
                SingleParameters param;

                param.epsilon = in_par.epsilon;
                param.sigma   = in_par.sigma;

                return param;

            }
        };

    }

    namespace Pairs{

        struct LennardJones {

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real epsilon;
                real sigma;
            };

            struct PairParameters{
                real epsilon;
                real sigma;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.epsilon = real(info.at("epsilon"));
                param.sigma   = real(info.at("sigma"));

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.epsilon = in_par.epsilon;
                param.sigma   = in_par.sigma;

                return param;
            }
        };

        struct Zhang {

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real radius;
                real epsilon;
                real mu;
                real chi;
                real theta;
            };

            struct PairParameters{
                real rb;
                real rc;
                real epsilon;
                real stheta;
                real mu;
                real chi;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.radius  = real(info.at("radius"));
                param.epsilon = real(info.at("epsilon"));
                param.mu      = real(info.at("mu"));
                param.chi     = real(info.at("chi"));
                param.theta   = real(info.at("theta"));

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.rb       = in_par.radius*pow(2,7.0/6.0);
                param.rc       = in_par.radius*5.2;
                param.epsilon  = in_par.epsilon;
                param.stheta   = sin(in_par.theta);
                param.mu       = in_par.mu;
                param.chi      = in_par.chi;

                return param;
            }
        };

        struct DistanceSwitchExponential {

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real E;
                real K;
                real rc;
            };

            struct PairParameters{
                real E;
                real K;
                real rc;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.E  = real(info.at("E"));
                param.K  = real(info.at("K"));
                param.rc = real(info.at("rc"));

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.E  = in_par.E;
                param.K  = in_par.K;
                param.rc = in_par.rc;

                return param;
            }
        };

        struct DistanceSwitchCosine {

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real E;
                real r_start;
                real rc;
            };

            struct PairParameters{
                real E;
                real r_start;
                real rc;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.E       = real(info.at("E"));
                param.r_start = real(info.at("r_start"));
                param.rc      = real(info.at("rc"));

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.E        = in_par.E;
                param.r_start  = in_par.r_start;
                param.rc       = in_par.rc;

                return param;
            }
        };

        template<typename potential>
        struct Helix {

            static inline __host__ tensor3 loadR_H(real3 e_x,
                                                   real3 e_y,
                                                   real3 e_z){

                // Check e_x, e_y, e_z are orthonormal
                real e_x_norm = sqrt(dot(e_x,e_x));
                real e_y_norm = sqrt(dot(e_y,e_y));
                real e_z_norm = sqrt(dot(e_z,e_z));

                if (abs(e_x_norm-1.0) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_x is not a unit vector");
                }
                if (abs(e_y_norm-1.0) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_y is not a unit vector");
                }
                if (abs(e_z_norm-1.0) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_z is not a unit vector");
                }

                real e_xy = dot(e_x,e_y);
                real e_xz = dot(e_x,e_z);
                real e_yz = dot(e_y,e_z);

                if (abs(e_xy) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_x and e_y are not orthogonal");
                }
                if (abs(e_xz) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_x and e_z are not orthogonal");
                }
                if (abs(e_yz) > 1e-6) {
                    System::log<System::CRITICAL>("[Helix] e_y and e_z are not orthogonal");
                }

                tensor3 R_H;

                R_H.xx = e_x.x;
                R_H.yx = e_x.y;
                R_H.zx = e_x.z;

                R_H.xy = e_y.x;
                R_H.yy = e_y.y;
                R_H.zy = e_y.z;

                R_H.xz= e_z.x;
                R_H.yz= e_z.y;
                R_H.zz= e_z.z;

                return R_H;
            }

            struct InputPairParameters{

                std::string name_i;
                std::string name_j;

                real Eb;
                typename potential::params helixParams;

                tensor3 R_H;
            };

            struct PairParameters{

                real Eb;
                typename potential::params helixParams;

                tensor3 R_H;
            };

            template<typename T>
            static inline __host__ InputPairParameters readPairParameters(std::map<std::string,T>& info){

                InputPairParameters param;

                param.name_i = std::string(info.at("name_i"));
                param.name_j = std::string(info.at("name_j"));

                param.Eb  = real(info.at("Eb"));

                param.helixParams = potential::readParamsMap(info);

                real3 e_x = real3(info.at("e_x"));
                real3 e_y = real3(info.at("e_y"));
                real3 e_z = real3(info.at("e_z"));

                param.R_H = loadR_H(e_x,e_y,e_z);

                System::log<System::MESSAGE>("[HelixParameters] Added parameter Eb: %f",param.Eb);

                System::log<System::MESSAGE>("              %f,%f,%f",param.R_H.xx,param.R_H.xy,param.R_H.xz);
                System::log<System::MESSAGE>("[Helix] R_H = %f,%f,%f",param.R_H.yx,param.R_H.yy,param.R_H.yz);
                System::log<System::MESSAGE>("              %f,%f,%f",param.R_H.zx,param.R_H.zy,param.R_H.zz);

                return param;

            }

            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){

                PairParameters param;

                param.Eb  = in_par.Eb;

                param.helixParams = in_par.helixParams;

                param.R_H = in_par.R_H;

                return param;
            }
        };

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

#endif
