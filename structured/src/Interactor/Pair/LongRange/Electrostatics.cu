#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"

#include "Interactor/Interactor.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Definitions/SFINAE.cuh"
#include "Utils/Containers/SetUtils.cuh"

#include "Definitions/Computations.cuh"
#include "Definitions/Types.cuh"

#include "Interactor/SpectralEwaldPoisson.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    class Electrostatics: public Interactor{

        protected:

            std::shared_ptr<GlobalData>     gd;

            ////////////////////////////////////

            std::unique_ptr<Poisson>        poisson;

            Box box;

            ////////////////////////////////////
            //Warnings

            bool warningEnergy           = false;
            bool warningForce            = false;
            bool warningLambdaDerivative = false;
            bool warningStress           = false;
            bool warningHessian          = false;
            bool warningMagneticField    = false;
            bool warningPairwiseForces   = false;

        public:

            Electrostatics(std::shared_ptr<GlobalData>           gd,
                           std::shared_ptr<ParticleGroup>        pg,
                           DataEntry& data,
                           std::string name):Interactor(pg,"Electrostatics: \"" +name+"\""),
                                             gd(gd){
                box = gd->getEnsemble()->getBox();

                // Read input parameters
                real dielectricConstant = data.getParameter<real>("dielectricConstant");
                real gw                 = data.getParameter<real>("gaussianWidth");
                real tolerance          = data.getParameter<real>("tolerance",1e-4);
                real splitThreshold     = data.getParameter<real>("splitThreshold",0.0);

                real eps = 1.0/gd->getUnits()->getElectricConversionFactor();
                // now eps = 4*pi*eps0 but eps = eps0*dielectricConstant
                eps *= dielectricConstant/(4.0*M_PI);

                Poisson::Parameters params;
                params.box = box;
                params.epsilon = eps;
                params.gw = gw;
                params.tolerance = tolerance;
                if (splitThreshold != 0.0){
                    params.split = splitThreshold;
                }

                poisson = std::make_unique<Poisson>(pg,params);
            }

            void sum(Computables comp,cudaStream_t st) override {

                Box box = gd->getEnsemble()->getBox();
                if (box != this->box){
                    System::log<System::CRITICAL>("[Electrostatics] The box has changed. This is not allowed.");
                }

                // Only energy and force computables are implemented
                if(comp.energy == true){
                    // This is done to force to compute only one computable at a time
                    Computables compEnergy;
                    compEnergy.energy = true;

                    poisson->sum(compEnergy,st);
                }

                if(comp.force == true){
                    // This is done to force to compute only one computable at a time
                    Computables compForce;
                    compForce.force = true;

                    poisson->sum(compForce,st);
                }

                if(comp.lambdaDerivative == true){
                    if(!warningLambdaDerivative){
                        System::log<System::WARNING>("[Electrostatics] (%s) Requested non-implemented computable (lambdaDerivative)",
                                                     name.c_str());
                        warningLambdaDerivative = true;
                    }
                }

                if(comp.stress == true){
                    if(!warningStress){
                        System::log<System::WARNING>("[Electrostatics] (%s) Requested non-implemented computable (stress)",
                                                     name.c_str());
                        warningStress = true;
                    }
                }

                if(comp.hessian == true){
                    if(!warningHessian){
                        System::log<System::WARNING>("[Electrostatics] (%s) Requested non-implemented computable (hessian)",
                                                     name.c_str());
                        warningHessian = true;
                    }
                }

                if(comp.magneticField == true){
                    if(!warningMagneticField){
                        System::log<System::WARNING>("[Electrostatics] (%s) Requested non-implemented computable (magneticField)",
                                                     name.c_str());
                        warningMagneticField = true;
                    }
                }

                if(comp.pairwiseForce == true){
                    if(!warningPairwiseForces){
                        System::log<System::WARNING>("[Electrostatics] (%s) Requested non-implemented computable (pairwiseForces)",
                                                     name.c_str());
                        warningPairwiseForces = true;
                    }
                }


            }
    };

}}}

REGISTER_INTERACTOR(
    LongRange,Electrostatics,
    uammd::structured::Interactor::Electrostatics
)
