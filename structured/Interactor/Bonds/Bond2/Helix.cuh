#ifndef __HELIX_BOND2__
#define __HELIX_BOND2__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    template<typename potential>
    struct Helix_{

        struct ComputationalData{
            real4* pos;
            real4* dir;

            Box    box;

            real Kb;
            real E;

            typename potential::params helixParams;

            tensor3 R_H;

            real3 e_next;
            real3 e_prev;

        };

        //Potential parameters

        struct StorageData{
            real Kb;
            typename potential::params helixParams;

            real E;

            tensor3 R_H;

            real3 e_next;
            real3 e_prev;
        };

        struct BondParameters{};

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                               std::shared_ptr<ParticleGroup>        pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.dir = pd->getDir(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            computational.Kb     = storage.Kb;
            computational.helixParams = storage.helixParams;

            computational.E = storage.E;

            computational.R_H = storage.R_H;

            computational.e_next = storage.e_next;
            computational.e_prev = storage.e_prev;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.Kb = data.getParameter<real>("Kb");
            System::log<System::MESSAGE>("[Helix] Kb = %f",storage.Kb);

            storage.helixParams = potential::readParams(data);

            storage.E = data.getParameter<real>("E");
            System::log<System::MESSAGE>("[Helix] E = %f",storage.E);

            real3 e_x = data.getParameter<real3>("e_x");
            real3 e_y = data.getParameter<real3>("e_y");
            real3 e_z = data.getParameter<real3>("e_z");

            storage.R_H = BasicParameters::Pairs::Helix<potential>::loadR_H(e_x,e_y,e_z);

            System::log<System::MESSAGE>("              %f,%f,%f",storage.R_H.xx,storage.R_H.xy,storage.R_H.xz);
            System::log<System::MESSAGE>("[Helix] R_H = %f,%f,%f",storage.R_H.yx,storage.R_H.yy,storage.R_H.yz);
            System::log<System::MESSAGE>("              %f,%f,%f",storage.R_H.zx,storage.R_H.zy,storage.R_H.zz);

            storage.e_next = data.getParameter<real3>("e_next");
            storage.e_prev = data.getParameter<real3>("e_prev");

            //Note that e_next and e_prev are not necessarily unit vectors

            System::log<System::MESSAGE>("[Helix] e_next = %f,%f,%f",storage.e_next.x,storage.e_next.y,storage.e_next.z);
            System::log<System::MESSAGE>("[Helix] e_prev = %f,%f,%f",storage.e_prev.x,storage.e_prev.y,storage.e_prev.z);

            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;
            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            int index_s;
            int index_e;

            real3 pp_s;
            real3 pp_e;

            tensor3 R_H = computational.R_H;

            if (currentParticleIndex == index_i){
                index_s = index_i;
                index_e = index_j;

                pp_s = computational.e_next;
                pp_e = computational.e_prev;
            } else {
                index_s = index_j;
                index_e = index_i;

                pp_s = computational.e_prev;
                pp_e = computational.e_next;

                R_H = R_H.transpose();
            }

            const Quat q_s = computational.dir[index_s];
            const Quat q_e = computational.dir[index_e];

            const real3 pos_s = make_real3(computational.pos[index_s]) + rotateVector(q_s,pp_s);
            const real3 pos_e = make_real3(computational.pos[index_e]) + rotateVector(q_e,pp_e);

            const real3 drse = computational.box.apply_pbc(pos_e - pos_s);

            const real Kb = computational.Kb;
            const typename potential::params p = computational.helixParams;

            const real E = computational.E;

            real e = BasicPotentials::Helix::Fixed::energy<potential>(drse,
                                                                      q_s,q_e,
                                                                      R_H,
                                                                      Kb,
                                                                      p.thetaParams,
                                                                      p.phiParams,
                                                                      E);

            return e;
        }

        static inline __device__ ForceTorque forceTorque(int index_i, int index_j,
                                                         int currentParticleIndex,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){

            int index_s;
            int index_e;

            real3 pp_s;
            real3 pp_e;

            tensor3 R_H = computational.R_H;

            if (currentParticleIndex == index_i){
                index_s = index_i;
                index_e = index_j;

                pp_s = computational.e_next;
                pp_e = computational.e_prev;
            } else {
                index_s = index_j;
                index_e = index_i;

                pp_s = computational.e_prev;
                pp_e = computational.e_next;

                R_H = R_H.transpose();
            }

            const Quat q_s = computational.dir[index_s];
            const Quat q_e = computational.dir[index_e];

            const real3 lps = rotateVector(q_s,pp_s);
            const real3 pos_s = make_real3(computational.pos[index_s]) + lps;
            const real3 pos_e = make_real3(computational.pos[index_e]) + rotateVector(q_e,pp_e);

            const real3 drse = computational.box.apply_pbc(pos_e - pos_s);

            const real Kb = computational.Kb;
            const typename potential::params p = computational.helixParams;

            const real E = computational.E;

            ForceTorque fT;

            fT.force  = make_real4(BasicPotentials::Helix::Fixed::force(drse,
                                                                        Kb),real(0.0));

            fT.torque = make_real4(BasicPotentials::Helix::Fixed::torque<potential>(drse,
                                                                                    lps,
                                                                                    q_s,q_e,
                                                                                    R_H,
                                                                                    Kb,
                                                                                    p.thetaParams,
                                                                                    p.phiParams,
                                                                                    E),real(0.0));

            return fT;
        }


    };

    using HelixExponential = Bond2Torque_<Helix_<typename BasicPotentials::Helix::exponential>>;
    using HelixCosine      = Bond2Torque_<Helix_<typename BasicPotentials::Helix::cosine>>;

}}}}

#endif
