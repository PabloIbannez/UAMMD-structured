#ifndef __GENERIC_LOADER__
#define __GENERIC_LOADER__
namespace uammd{
namespace structured{
namespace Potentials{
namespace GenericLoader{

    bool isInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();
        
        if("Bond1" == potType and "FixedHarmonic" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Fene" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "FeneCommon_K_R0" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "FeneCommon_r0_K_R0" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "DebyeHuckel" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Harmonic" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LambdaHarmonic" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "HarmonicCommon_K" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "HarmonicCommon_K_r0" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Steric6" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Steric6Common_epsilon_sigma" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Steric12" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Steric12Common_epsilon_sigma" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Morse" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "MorseCommon_D" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "MorseCommon_r0_E_D" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "MorseWCACommon_eps0" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesType1" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesType1Common_epsilon" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesType2" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesType2Common_epsilon" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesType3" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesType3Common_epsilon" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesKaranicolasBrooks" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesKaranicolasBrooksCommon_epsilon" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesGaussian" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesGaussianCommon_epsilon_D" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesSoftCoreType1" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesSoftCoreType2" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesSoftCoreType1Common_epsilon" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "LennardJonesSoftCoreType2Common_epsilon" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "Gaussian" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "GaussianCommon_E_r0_D" == potSubType){
            return true;
        }
            
        if("Bond3" == potType and "BestChenHummerAngular" == potSubType){
            return true;
        }
            
        if("Bond3" == potType and "KratkyPorod" == potSubType){
            return true;
        }
            
        if("Bond3" == potType and "KratkyPorodCommon_K" == potSubType){
            return true;
        }
            
        if("Bond3" == potType and "HarmonicAngular" == potSubType){
            return true;
        }
            
        if("Bond3" == potType and "HarmonicAngularCommon_K" == potSubType){
            return true;
        }
            
        if("Bond3" == potType and "HarmonicAngularCommon_K_theta0" == potSubType){
            return true;
        }
            
        if("Bond4" == potType and "Dihedral" == potSubType){
            return true;
        }
            
        if("Bond4" == potType and "DihedralCommon_n_K_phi0" == potSubType){
            return true;
        }
            
        if("Bond4" == potType and "Dihedral4" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "MaxDistanceRestraint" == potSubType){
            return true;
        }
            
        if("Bond1" == potType and "ConstantForce" == potSubType){
            return true;
        }
            
        if("Bond1" == potType and "LambdaFixedHarmonicAnisotropic" == potSubType){
            return true;
        }
            
        if("Bond4" == potType and "IDP_Fourier" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "HelixExponential" == potSubType){
            return true;
        }
            
        if("Bond2" == potType and "HelixCosine" == potSubType){
            return true;
        }
            
        if("Bond1" == potType and "FixedHarmonicAnisotropic" == potSubType){
            return true;
        }
            
        if("External" == potType and "ConstantForce" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "SurfaceLennardJonesType1" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "SurfaceLennardJonesType2" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "SurfaceWCAType1" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "SurfaceWCAType2" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "SurfaceGeneralLennardJonesType1" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "SurfaceGeneralLennardJonesType2" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "SurfaceAnchorage" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "ParabolaSurface" == potSubType){
            return true;
        }
            
        if("External" == potType and "ACMagneticField" == potSubType){
            return true;
        }
            
        if("External" == potType and "ConstantMagneticField" == potSubType){
            return true;
        }
            
        if("External" == potType and "SphericalShell" == potSubType){
            return true;
        }
            
        if("Surface" == potType and "Absorbed" == potSubType){
            return true;
        }
            
        if("External" == potType and "Plates" == potSubType){
            return true;
        }
            
        if("External" == potType and "ConstantTorque" == potSubType){
            return true;
        }
            
        if("Set1" == potType and "ConstantTorqueOverCenterOfMass" == potSubType){
            return true;
        }
            
        if("Set1" == potType and "ConstantForceOverCenterOfMass" == potSubType){
            return true;
        }
            
        if("Set1" == potType and "FixedHarmonicCenterOfMass" == potSubType){
            return true;
        }
            
        if("Set1" == potType and "FixedHarmonicCommon_K_r0CenterOfMass" == potSubType){
            return true;
        }
            
        if("Set2" == potType and "ConstantTorqueBetweenCentersOfMass" == potSubType){
            return true;
        }
            
        if("Set2" == potType and "ConstantForceBetweenCentersOfMass" == potSubType){
            return true;
        }
            
        if("Set2" == potType and "HarmonicBondBetweenCentersOfMass" == potSubType){
            return true;
        }
            
        if("AFM" == potType and "SphericalTip" == potSubType){
            return true;
        }
            
        if("AFM" == potType and "SphericallyBluntedConicTip" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "DebyeHuckel" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "Zhang" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "DLVOType1" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "DLVOType2" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "DLVOType3" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "Clashed" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "KimHummer" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "LennardJonesType1" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "LennardJonesType2" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "LennardJonesType3" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "WCAType1" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "WCAType2" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "WCAType3" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "GeneralLennardJonesType1" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "GeneralLennardJonesType2" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "GeneralLennardJonesType3" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "Steric6" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "Steric12" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "StericInner12" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "StericInner6" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "DipolarMagnetic" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "Steric6SoftCore" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "Steric12SoftCore" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "LennardJonesSoftCoreType1" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "LennardJonesSoftCoreType2" == potSubType){
            return true;
        }
            
        if("NonBonded" == potType and "SplitLennardJones" == potSubType){
            return true;
        }
            
        if("PatchyParticles" == potType and "PatchyParticles" == potSubType){
            return true;
        }
            
        if("PatchyParticles" == potType and "DynamicallyBondedPatchyParticles" == potSubType){
            return true;
        }
            
        return false;

    }

    
    std::shared_ptr<typename uammd::Interactor>
    loadGeneric(std::shared_ptr<ExtendedSystem> sys,
                std::shared_ptr<GlobalData>     gd,
                std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        std::shared_ptr<typename uammd::Interactor> interactor;
        bool found = false;

        if("Bond1" == potType and "FixedHarmonic" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond1::FixedHarmonic potential",path.back().c_str());

            std::shared_ptr<Bond1::FixedHarmonic> pot =
            std::make_shared<Bond1::FixedHarmonic>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond1::FixedHarmonic>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Fene" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Fene potential",path.back().c_str());

            std::shared_ptr<Bond2::Fene> pot =
            std::make_shared<Bond2::Fene>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Fene>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "FeneCommon_K_R0" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::FeneCommon_K_R0 potential",path.back().c_str());

            std::shared_ptr<Bond2::FeneCommon_K_R0> pot =
            std::make_shared<Bond2::FeneCommon_K_R0>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::FeneCommon_K_R0>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "FeneCommon_r0_K_R0" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::FeneCommon_r0_K_R0 potential",path.back().c_str());

            std::shared_ptr<Bond2::FeneCommon_r0_K_R0> pot =
            std::make_shared<Bond2::FeneCommon_r0_K_R0>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::FeneCommon_r0_K_R0>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "DebyeHuckel" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::DebyeHuckel potential",path.back().c_str());

            std::shared_ptr<Bond2::DebyeHuckel> pot =
            std::make_shared<Bond2::DebyeHuckel>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::DebyeHuckel>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Harmonic" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Harmonic potential",path.back().c_str());

            std::shared_ptr<Bond2::Harmonic> pot =
            std::make_shared<Bond2::Harmonic>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Harmonic>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LambdaHarmonic" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LambdaHarmonic potential",path.back().c_str());

            std::shared_ptr<Bond2::LambdaHarmonic> pot =
            std::make_shared<Bond2::LambdaHarmonic>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LambdaHarmonic>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "HarmonicCommon_K" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::HarmonicCommon_K potential",path.back().c_str());

            std::shared_ptr<Bond2::HarmonicCommon_K> pot =
            std::make_shared<Bond2::HarmonicCommon_K>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::HarmonicCommon_K>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "HarmonicCommon_K_r0" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::HarmonicCommon_K_r0 potential",path.back().c_str());

            std::shared_ptr<Bond2::HarmonicCommon_K_r0> pot =
            std::make_shared<Bond2::HarmonicCommon_K_r0>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::HarmonicCommon_K_r0>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Steric6" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Steric6 potential",path.back().c_str());

            std::shared_ptr<Bond2::Steric6> pot =
            std::make_shared<Bond2::Steric6>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Steric6>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Steric6Common_epsilon_sigma" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Steric6Common_epsilon_sigma potential",path.back().c_str());

            std::shared_ptr<Bond2::Steric6Common_epsilon_sigma> pot =
            std::make_shared<Bond2::Steric6Common_epsilon_sigma>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Steric6Common_epsilon_sigma>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Steric12" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Steric12 potential",path.back().c_str());

            std::shared_ptr<Bond2::Steric12> pot =
            std::make_shared<Bond2::Steric12>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Steric12>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Steric12Common_epsilon_sigma" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Steric12Common_epsilon_sigma potential",path.back().c_str());

            std::shared_ptr<Bond2::Steric12Common_epsilon_sigma> pot =
            std::make_shared<Bond2::Steric12Common_epsilon_sigma>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Steric12Common_epsilon_sigma>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Morse" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Morse potential",path.back().c_str());

            std::shared_ptr<Bond2::Morse> pot =
            std::make_shared<Bond2::Morse>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Morse>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "MorseCommon_D" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::MorseCommon_D potential",path.back().c_str());

            std::shared_ptr<Bond2::MorseCommon_D> pot =
            std::make_shared<Bond2::MorseCommon_D>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::MorseCommon_D>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "MorseCommon_r0_E_D" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::MorseCommon_r0_E_D potential",path.back().c_str());

            std::shared_ptr<Bond2::MorseCommon_r0_E_D> pot =
            std::make_shared<Bond2::MorseCommon_r0_E_D>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::MorseCommon_r0_E_D>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "MorseWCACommon_eps0" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::MorseWCACommon_eps0 potential",path.back().c_str());

            std::shared_ptr<Bond2::MorseWCACommon_eps0> pot =
            std::make_shared<Bond2::MorseWCACommon_eps0>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::MorseWCACommon_eps0>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesType1 potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesType1> pot =
            std::make_shared<Bond2::LennardJonesType1>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesType1>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesType1Common_epsilon" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesType1Common_epsilon potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesType1Common_epsilon> pot =
            std::make_shared<Bond2::LennardJonesType1Common_epsilon>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesType1Common_epsilon>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesType2 potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesType2> pot =
            std::make_shared<Bond2::LennardJonesType2>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesType2>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesType2Common_epsilon" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesType2Common_epsilon potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesType2Common_epsilon> pot =
            std::make_shared<Bond2::LennardJonesType2Common_epsilon>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesType2Common_epsilon>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesType3" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesType3 potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesType3> pot =
            std::make_shared<Bond2::LennardJonesType3>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesType3>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesType3Common_epsilon" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesType3Common_epsilon potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesType3Common_epsilon> pot =
            std::make_shared<Bond2::LennardJonesType3Common_epsilon>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesType3Common_epsilon>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesKaranicolasBrooks" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesKaranicolasBrooks potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesKaranicolasBrooks> pot =
            std::make_shared<Bond2::LennardJonesKaranicolasBrooks>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesKaranicolasBrooks>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesKaranicolasBrooksCommon_epsilon" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesKaranicolasBrooksCommon_epsilon potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesKaranicolasBrooksCommon_epsilon> pot =
            std::make_shared<Bond2::LennardJonesKaranicolasBrooksCommon_epsilon>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesKaranicolasBrooksCommon_epsilon>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesGaussian" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesGaussian potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesGaussian> pot =
            std::make_shared<Bond2::LennardJonesGaussian>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesGaussian>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesGaussianCommon_epsilon_D" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesGaussianCommon_epsilon_D potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesGaussianCommon_epsilon_D> pot =
            std::make_shared<Bond2::LennardJonesGaussianCommon_epsilon_D>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesGaussianCommon_epsilon_D>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesSoftCoreType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesSoftCoreType1 potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesSoftCoreType1> pot =
            std::make_shared<Bond2::LennardJonesSoftCoreType1>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesSoftCoreType1>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesSoftCoreType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesSoftCoreType2 potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesSoftCoreType2> pot =
            std::make_shared<Bond2::LennardJonesSoftCoreType2>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesSoftCoreType2>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesSoftCoreType1Common_epsilon" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesSoftCoreType1Common_epsilon potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesSoftCoreType1Common_epsilon> pot =
            std::make_shared<Bond2::LennardJonesSoftCoreType1Common_epsilon>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesSoftCoreType1Common_epsilon>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "LennardJonesSoftCoreType2Common_epsilon" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::LennardJonesSoftCoreType2Common_epsilon potential",path.back().c_str());

            std::shared_ptr<Bond2::LennardJonesSoftCoreType2Common_epsilon> pot =
            std::make_shared<Bond2::LennardJonesSoftCoreType2Common_epsilon>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::LennardJonesSoftCoreType2Common_epsilon>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "Gaussian" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::Gaussian potential",path.back().c_str());

            std::shared_ptr<Bond2::Gaussian> pot =
            std::make_shared<Bond2::Gaussian>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::Gaussian>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "GaussianCommon_E_r0_D" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::GaussianCommon_E_r0_D potential",path.back().c_str());

            std::shared_ptr<Bond2::GaussianCommon_E_r0_D> pot =
            std::make_shared<Bond2::GaussianCommon_E_r0_D>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::GaussianCommon_E_r0_D>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond3" == potType and "BestChenHummerAngular" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond3::BestChenHummerAngular potential",path.back().c_str());

            std::shared_ptr<Bond3::BestChenHummerAngular> pot =
            std::make_shared<Bond3::BestChenHummerAngular>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond3::BestChenHummerAngular>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond3" == potType and "KratkyPorod" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond3::KratkyPorod potential",path.back().c_str());

            std::shared_ptr<Bond3::KratkyPorod> pot =
            std::make_shared<Bond3::KratkyPorod>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond3::KratkyPorod>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond3" == potType and "KratkyPorodCommon_K" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond3::KratkyPorodCommon_K potential",path.back().c_str());

            std::shared_ptr<Bond3::KratkyPorodCommon_K> pot =
            std::make_shared<Bond3::KratkyPorodCommon_K>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond3::KratkyPorodCommon_K>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond3" == potType and "HarmonicAngular" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond3::HarmonicAngular potential",path.back().c_str());

            std::shared_ptr<Bond3::HarmonicAngular> pot =
            std::make_shared<Bond3::HarmonicAngular>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond3::HarmonicAngular>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond3" == potType and "HarmonicAngularCommon_K" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond3::HarmonicAngularCommon_K potential",path.back().c_str());

            std::shared_ptr<Bond3::HarmonicAngularCommon_K> pot =
            std::make_shared<Bond3::HarmonicAngularCommon_K>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond3::HarmonicAngularCommon_K>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond3" == potType and "HarmonicAngularCommon_K_theta0" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond3::HarmonicAngularCommon_K_theta0 potential",path.back().c_str());

            std::shared_ptr<Bond3::HarmonicAngularCommon_K_theta0> pot =
            std::make_shared<Bond3::HarmonicAngularCommon_K_theta0>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond3::HarmonicAngularCommon_K_theta0>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond4" == potType and "Dihedral" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond4::Dihedral potential",path.back().c_str());

            std::shared_ptr<Bond4::Dihedral> pot =
            std::make_shared<Bond4::Dihedral>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond4::Dihedral>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond4" == potType and "DihedralCommon_n_K_phi0" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond4::DihedralCommon_n_K_phi0 potential",path.back().c_str());

            std::shared_ptr<Bond4::DihedralCommon_n_K_phi0> pot =
            std::make_shared<Bond4::DihedralCommon_n_K_phi0>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond4::DihedralCommon_n_K_phi0>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond4" == potType and "Dihedral4" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond4::Dihedral4 potential",path.back().c_str());

            std::shared_ptr<Bond4::Dihedral4> pot =
            std::make_shared<Bond4::Dihedral4>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond4::Dihedral4>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "MaxDistanceRestraint" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::MaxDistanceRestraint potential",path.back().c_str());

            std::shared_ptr<Bond2::MaxDistanceRestraint> pot =
            std::make_shared<Bond2::MaxDistanceRestraint>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::MaxDistanceRestraint>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond1" == potType and "ConstantForce" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond1::ConstantForce potential",path.back().c_str());

            std::shared_ptr<Bond1::ConstantForce> pot =
            std::make_shared<Bond1::ConstantForce>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond1::ConstantForce>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond1" == potType and "LambdaFixedHarmonicAnisotropic" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond1::LambdaFixedHarmonicAnisotropic potential",path.back().c_str());

            std::shared_ptr<Bond1::LambdaFixedHarmonicAnisotropic> pot =
            std::make_shared<Bond1::LambdaFixedHarmonicAnisotropic>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond1::LambdaFixedHarmonicAnisotropic>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond4" == potType and "IDP_Fourier" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond4::IDP_Fourier potential",path.back().c_str());

            std::shared_ptr<Bond4::IDP_Fourier> pot =
            std::make_shared<Bond4::IDP_Fourier>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond4::IDP_Fourier>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "HelixExponential" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::HelixExponential potential",path.back().c_str());

            std::shared_ptr<Bond2::HelixExponential> pot =
            std::make_shared<Bond2::HelixExponential>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::HelixExponential>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond2" == potType and "HelixCosine" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond2::HelixCosine potential",path.back().c_str());

            std::shared_ptr<Bond2::HelixCosine> pot =
            std::make_shared<Bond2::HelixCosine>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond2::HelixCosine>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Bond1" == potType and "FixedHarmonicAnisotropic" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Bond1::FixedHarmonicAnisotropic potential",path.back().c_str());

            std::shared_ptr<Bond1::FixedHarmonicAnisotropic> pot =
            std::make_shared<Bond1::FixedHarmonicAnisotropic>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::BondsInteractor<Bond1::FixedHarmonicAnisotropic>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("External" == potType and "ConstantForce" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected External::ConstantForce potential",path.back().c_str());

            std::shared_ptr<External::ConstantForce> pot =
            std::make_shared<External::ConstantForce>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<External::ConstantForce>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "SurfaceLennardJonesType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::SurfaceLennardJonesType1 potential",path.back().c_str());

            std::shared_ptr<Surface::SurfaceLennardJonesType1> pot =
            std::make_shared<Surface::SurfaceLennardJonesType1>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::SurfaceLennardJonesType1>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "SurfaceLennardJonesType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::SurfaceLennardJonesType2 potential",path.back().c_str());

            std::shared_ptr<Surface::SurfaceLennardJonesType2> pot =
            std::make_shared<Surface::SurfaceLennardJonesType2>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::SurfaceLennardJonesType2>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "SurfaceWCAType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::SurfaceWCAType1 potential",path.back().c_str());

            std::shared_ptr<Surface::SurfaceWCAType1> pot =
            std::make_shared<Surface::SurfaceWCAType1>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::SurfaceWCAType1>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "SurfaceWCAType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::SurfaceWCAType2 potential",path.back().c_str());

            std::shared_ptr<Surface::SurfaceWCAType2> pot =
            std::make_shared<Surface::SurfaceWCAType2>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::SurfaceWCAType2>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "SurfaceGeneralLennardJonesType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::SurfaceGeneralLennardJonesType1 potential",path.back().c_str());

            std::shared_ptr<Surface::SurfaceGeneralLennardJonesType1> pot =
            std::make_shared<Surface::SurfaceGeneralLennardJonesType1>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::SurfaceGeneralLennardJonesType1>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "SurfaceGeneralLennardJonesType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::SurfaceGeneralLennardJonesType2 potential",path.back().c_str());

            std::shared_ptr<Surface::SurfaceGeneralLennardJonesType2> pot =
            std::make_shared<Surface::SurfaceGeneralLennardJonesType2>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::SurfaceGeneralLennardJonesType2>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "SurfaceAnchorage" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::SurfaceAnchorage potential",path.back().c_str());

            std::shared_ptr<Surface::SurfaceAnchorage> pot =
            std::make_shared<Surface::SurfaceAnchorage>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::SurfaceAnchorage>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "ParabolaSurface" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::ParabolaSurface potential",path.back().c_str());

            std::shared_ptr<Surface::ParabolaSurface> pot =
            std::make_shared<Surface::ParabolaSurface>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::ParabolaSurface>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("External" == potType and "ACMagneticField" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected External::ACMagneticField potential",path.back().c_str());

            std::shared_ptr<External::ACMagneticField> pot =
            std::make_shared<External::ACMagneticField>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<External::ACMagneticField>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("External" == potType and "ConstantMagneticField" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected External::ConstantMagneticField potential",path.back().c_str());

            std::shared_ptr<External::ConstantMagneticField> pot =
            std::make_shared<External::ConstantMagneticField>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<External::ConstantMagneticField>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("External" == potType and "SphericalShell" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected External::SphericalShell potential",path.back().c_str());

            std::shared_ptr<External::SphericalShell> pot =
            std::make_shared<External::SphericalShell>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<External::SphericalShell>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Surface" == potType and "Absorbed" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Surface::Absorbed potential",path.back().c_str());

            std::shared_ptr<Surface::Absorbed> pot =
            std::make_shared<Surface::Absorbed>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<Surface::Absorbed>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("External" == potType and "Plates" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected External::Plates potential",path.back().c_str());

            std::shared_ptr<External::Plates> pot =
            std::make_shared<External::Plates>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<External::Plates>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("External" == potType and "ConstantTorque" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected External::ConstantTorque potential",path.back().c_str());

            std::shared_ptr<External::ConstantTorque> pot =
            std::make_shared<External::ConstantTorque>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SingleInteractor<External::ConstantTorque>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Set1" == potType and "ConstantTorqueOverCenterOfMass" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Set1::ConstantTorqueOverCenterOfMass potential",path.back().c_str());

            std::shared_ptr<Set1::ConstantTorqueOverCenterOfMass> pot =
            std::make_shared<Set1::ConstantTorqueOverCenterOfMass>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SetInteractor<Set1::ConstantTorqueOverCenterOfMass>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Set1" == potType and "ConstantForceOverCenterOfMass" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Set1::ConstantForceOverCenterOfMass potential",path.back().c_str());

            std::shared_ptr<Set1::ConstantForceOverCenterOfMass> pot =
            std::make_shared<Set1::ConstantForceOverCenterOfMass>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SetInteractor<Set1::ConstantForceOverCenterOfMass>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Set1" == potType and "FixedHarmonicCenterOfMass" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Set1::FixedHarmonicCenterOfMass potential",path.back().c_str());

            std::shared_ptr<Set1::FixedHarmonicCenterOfMass> pot =
            std::make_shared<Set1::FixedHarmonicCenterOfMass>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SetInteractor<Set1::FixedHarmonicCenterOfMass>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Set1" == potType and "FixedHarmonicCommon_K_r0CenterOfMass" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Set1::FixedHarmonicCommon_K_r0CenterOfMass potential",path.back().c_str());

            std::shared_ptr<Set1::FixedHarmonicCommon_K_r0CenterOfMass> pot =
            std::make_shared<Set1::FixedHarmonicCommon_K_r0CenterOfMass>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SetInteractor<Set1::FixedHarmonicCommon_K_r0CenterOfMass>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Set2" == potType and "ConstantTorqueBetweenCentersOfMass" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Set2::ConstantTorqueBetweenCentersOfMass potential",path.back().c_str());

            std::shared_ptr<Set2::ConstantTorqueBetweenCentersOfMass> pot =
            std::make_shared<Set2::ConstantTorqueBetweenCentersOfMass>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SetInteractor<Set2::ConstantTorqueBetweenCentersOfMass>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Set2" == potType and "ConstantForceBetweenCentersOfMass" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Set2::ConstantForceBetweenCentersOfMass potential",path.back().c_str());

            std::shared_ptr<Set2::ConstantForceBetweenCentersOfMass> pot =
            std::make_shared<Set2::ConstantForceBetweenCentersOfMass>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SetInteractor<Set2::ConstantForceBetweenCentersOfMass>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("Set2" == potType and "HarmonicBondBetweenCentersOfMass" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected Set2::HarmonicBondBetweenCentersOfMass potential",path.back().c_str());

            std::shared_ptr<Set2::HarmonicBondBetweenCentersOfMass> pot =
            std::make_shared<Set2::HarmonicBondBetweenCentersOfMass>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::SetInteractor<Set2::HarmonicBondBetweenCentersOfMass>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("AFM" == potType and "SphericalTip" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected AFM::SphericalTip potential",path.back().c_str());

            std::shared_ptr<AFM::SphericalTip> pot =
            std::make_shared<AFM::SphericalTip>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::AFMInteractor<AFM::SphericalTip>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("AFM" == potType and "SphericallyBluntedConicTip" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected AFM::SphericallyBluntedConicTip potential",path.back().c_str());

            std::shared_ptr<AFM::SphericallyBluntedConicTip> pot =
            std::make_shared<AFM::SphericallyBluntedConicTip>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::AFMInteractor<AFM::SphericallyBluntedConicTip>>(gd,pg,data,pot,path.back());
            found = true;
        }
        if("NonBonded" == potType and "DebyeHuckel" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::DebyeHuckel potential",path.back().c_str());

            std::shared_ptr<NonBonded::DebyeHuckel> pot =
            std::make_shared<NonBonded::DebyeHuckel>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::DebyeHuckel,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "Zhang" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::Zhang potential",path.back().c_str());

            std::shared_ptr<NonBonded::Zhang> pot =
            std::make_shared<NonBonded::Zhang>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::Zhang,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "DLVOType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::DLVOType1 potential",path.back().c_str());

            std::shared_ptr<NonBonded::DLVOType1> pot =
            std::make_shared<NonBonded::DLVOType1>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::DLVOType1,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "DLVOType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::DLVOType2 potential",path.back().c_str());

            std::shared_ptr<NonBonded::DLVOType2> pot =
            std::make_shared<NonBonded::DLVOType2>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::DLVOType2,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "DLVOType3" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::DLVOType3 potential",path.back().c_str());

            std::shared_ptr<NonBonded::DLVOType3> pot =
            std::make_shared<NonBonded::DLVOType3>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::DLVOType3,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "Clashed" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::Clashed potential",path.back().c_str());

            std::shared_ptr<NonBonded::Clashed> pot =
            std::make_shared<NonBonded::Clashed>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::Clashed,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "KimHummer" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::KimHummer potential",path.back().c_str());

            std::shared_ptr<NonBonded::KimHummer> pot =
            std::make_shared<NonBonded::KimHummer>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::KimHummer,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "LennardJonesType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::LennardJonesType1 potential",path.back().c_str());

            std::shared_ptr<NonBonded::LennardJonesType1> pot =
            std::make_shared<NonBonded::LennardJonesType1>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::LennardJonesType1,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "LennardJonesType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::LennardJonesType2 potential",path.back().c_str());

            std::shared_ptr<NonBonded::LennardJonesType2> pot =
            std::make_shared<NonBonded::LennardJonesType2>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::LennardJonesType2,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "LennardJonesType3" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::LennardJonesType3 potential",path.back().c_str());

            std::shared_ptr<NonBonded::LennardJonesType3> pot =
            std::make_shared<NonBonded::LennardJonesType3>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::LennardJonesType3,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "WCAType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::WCAType1 potential",path.back().c_str());

            std::shared_ptr<NonBonded::WCAType1> pot =
            std::make_shared<NonBonded::WCAType1>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::WCAType1,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "WCAType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::WCAType2 potential",path.back().c_str());

            std::shared_ptr<NonBonded::WCAType2> pot =
            std::make_shared<NonBonded::WCAType2>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::WCAType2,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "WCAType3" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::WCAType3 potential",path.back().c_str());

            std::shared_ptr<NonBonded::WCAType3> pot =
            std::make_shared<NonBonded::WCAType3>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::WCAType3,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "GeneralLennardJonesType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::GeneralLennardJonesType1 potential",path.back().c_str());

            std::shared_ptr<NonBonded::GeneralLennardJonesType1> pot =
            std::make_shared<NonBonded::GeneralLennardJonesType1>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::GeneralLennardJonesType1,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "GeneralLennardJonesType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::GeneralLennardJonesType2 potential",path.back().c_str());

            std::shared_ptr<NonBonded::GeneralLennardJonesType2> pot =
            std::make_shared<NonBonded::GeneralLennardJonesType2>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::GeneralLennardJonesType2,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "GeneralLennardJonesType3" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::GeneralLennardJonesType3 potential",path.back().c_str());

            std::shared_ptr<NonBonded::GeneralLennardJonesType3> pot =
            std::make_shared<NonBonded::GeneralLennardJonesType3>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::GeneralLennardJonesType3,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "Steric6" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::Steric6 potential",path.back().c_str());

            std::shared_ptr<NonBonded::Steric6> pot =
            std::make_shared<NonBonded::Steric6>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::Steric6,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "Steric12" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::Steric12 potential",path.back().c_str());

            std::shared_ptr<NonBonded::Steric12> pot =
            std::make_shared<NonBonded::Steric12>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::Steric12,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "StericInner12" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::StericInner12 potential",path.back().c_str());

            std::shared_ptr<NonBonded::StericInner12> pot =
            std::make_shared<NonBonded::StericInner12>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::StericInner12,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "StericInner6" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::StericInner6 potential",path.back().c_str());

            std::shared_ptr<NonBonded::StericInner6> pot =
            std::make_shared<NonBonded::StericInner6>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::StericInner6,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "DipolarMagnetic" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::DipolarMagnetic potential",path.back().c_str());

            std::shared_ptr<NonBonded::DipolarMagnetic> pot =
            std::make_shared<NonBonded::DipolarMagnetic>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::DipolarMagnetic,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "Steric6SoftCore" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::Steric6SoftCore potential",path.back().c_str());

            std::shared_ptr<NonBonded::Steric6SoftCore> pot =
            std::make_shared<NonBonded::Steric6SoftCore>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::Steric6SoftCore,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "Steric12SoftCore" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::Steric12SoftCore potential",path.back().c_str());

            std::shared_ptr<NonBonded::Steric12SoftCore> pot =
            std::make_shared<NonBonded::Steric12SoftCore>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::Steric12SoftCore,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "LennardJonesSoftCoreType1" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::LennardJonesSoftCoreType1 potential",path.back().c_str());

            std::shared_ptr<NonBonded::LennardJonesSoftCoreType1> pot =
            std::make_shared<NonBonded::LennardJonesSoftCoreType1>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::LennardJonesSoftCoreType1,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "LennardJonesSoftCoreType2" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::LennardJonesSoftCoreType2 potential",path.back().c_str());

            std::shared_ptr<NonBonded::LennardJonesSoftCoreType2> pot =
            std::make_shared<NonBonded::LennardJonesSoftCoreType2>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::LennardJonesSoftCoreType2,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("NonBonded" == potType and "SplitLennardJones" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected NonBonded::SplitLennardJones potential",path.back().c_str());

            std::shared_ptr<NonBonded::SplitLennardJones> pot =
            std::make_shared<NonBonded::SplitLennardJones>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::PairInteractor<NonBonded::SplitLennardJones,VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }
        if("PatchyParticles" == potType and "PatchyParticles" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected PatchyParticles::PatchyParticles potential",path.back().c_str());
            interactor = std::make_shared<typename Interactor::PatchyParticles::PatchyParticles>(gd,pg,path,path.back());
            found = true;
        }
        if("PatchyParticles" == potType and "DynamicallyBondedPatchyParticles" == potSubType){
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected PatchyParticles::DynamicallyBondedPatchyParticles potential",path.back().c_str());
            interactor = std::make_shared<typename Interactor::PatchyParticles::DynamicallyBondedPatchyParticles>(gd,pg,path,path.back());
            found = true;
        }

        if(not found){
            System::log<System::CRITICAL>("[GenericLoader] (%s) Potential of class %s and subclass %s not found",
                                           path.back().c_str(),data.getType().c_str(),data.getSubType().c_str());
        }

        return interactor;

    }

    }}}}
#endif
