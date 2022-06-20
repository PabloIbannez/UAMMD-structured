
#ifndef __GENERIC__
#define __GENERIC__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Generic{

template<class Units_,
         class Types_,
         template <class Topology_> class Condition_>
class Generic : public ForceFieldBase<Units_,Types_>{
        
    protected:

        using Base = ForceFieldBase<Units_,Types_>;

        using Condition = Condition_<typename Base::Topology>;
        
        using NeighbourList = ConditionedVerletListSet<Condition>;
        
        std::shared_ptr<Condition> condition;
        std::shared_ptr<NeighbourList>    nl; 
            
        real VerletListDst;

        using Bond1HarmonicType = Potentials::Bond1::Harmonic;
        using Bond1HarmonicConst_r0Type = Potentials::Bond1::HarmonicConst_r0;
        using Bond1HarmonicConst_K_r0Type = Potentials::Bond1::HarmonicConst_K_r0;
        using Bond1HarmonicCommon_K_r0Type = Potentials::Bond1::HarmonicCommon_K_r0;
        using Bond1HarmonicConstCommon_K_r0Type = Potentials::Bond1::HarmonicConstCommon_K_r0;
        using Bond2FeneType = Potentials::Bond2::Fene;
        using Bond2FeneConst_K_R0Type = Potentials::Bond2::FeneConst_K_R0;
        using Bond2FeneConst_r0_K_R0Type = Potentials::Bond2::FeneConst_r0_K_R0;
        using Bond2DebyeHuckelType = Potentials::Bond2::DebyeHuckel<typename Base::Units>;
        using Bond2HarmonicType = Potentials::Bond2::Harmonic;
        using Bond2HarmonicConst_KType = Potentials::Bond2::HarmonicConst_K;
        using Bond2HarmonicConst_K_r0Type = Potentials::Bond2::HarmonicConst_K_r0;
        using Bond2Steric6Type = Potentials::Bond2::Steric6;
        using Bond2Steric6Const_epsilon_sigmaType = Potentials::Bond2::Steric6Const_epsilon_sigma;
        using Bond2Steric12Type = Potentials::Bond2::Steric12;
        using Bond2Steric12Const_epsilon_sigmaType = Potentials::Bond2::Steric12Const_epsilon_sigma;
        using Bond2MorseType = Potentials::Bond2::Morse;
        using Bond2MorseConst_DType = Potentials::Bond2::MorseConst_D;
        using Bond2MorseConst_r0_E_DType = Potentials::Bond2::MorseConst_r0_E_D;
        using Bond2MorseWCAType = Potentials::Bond2::MorseWCA;
        using Bond2LennardJonesType2Type = Potentials::Bond2::LennardJonesType2;
        using Bond2LennardJonesType2Const_eType = Potentials::Bond2::LennardJonesType2Const_e;
        using Bond2LennardJonesType3Type = Potentials::Bond2::LennardJonesType3;
        using Bond2LennardJonesType3Const_eType = Potentials::Bond2::LennardJonesType3Const_e;
        using Bond2LennardJonesKaranicolasBrooksType = Potentials::Bond2::LennardJonesKaranicolasBrooks;
        using Bond2LennardJonesGaussianType = Potentials::Bond2::LennardJonesGaussian;
        using Bond2LennardJonesGaussianConst_e_DType = Potentials::Bond2::LennardJonesGaussianConst_e_D;
        using Bond2GaussianType = Potentials::Bond2::Gaussian;
        using Bond2GaussianConst_E_r0_DType = Potentials::Bond2::GaussianConst_E_r0_D;
        using Bond2OrientedHarmonicType = Potentials::Bond2::OrientedHarmonic;
        using Bond3BestChenHummerAngularType = Potentials::Bond3::BestChenHummerAngular;
        using Bond3KratkyPorodType = Potentials::Bond3::KratkyPorod;
        using Bond3KratkyPorodConst_KType = Potentials::Bond3::KratkyPorodConst_K;
        using Bond3HarmonicAngularType = Potentials::Bond3::HarmonicAngular;
        using Bond3HarmonicAngularConst_KType = Potentials::Bond3::HarmonicAngularConst_K;
        using Bond3HarmonicAngularConst_K_ang0Type = Potentials::Bond3::HarmonicAngularConst_K_ang0;
        using Bond4Dihedral4Type = Potentials::Bond4::Dihedral4;
        using Bond4DihedralType = Potentials::Bond4::Dihedral;
        using Bond4DihedralConst_n_K_phi0Type = Potentials::Bond4::DihedralConst_n_K_phi0;
        using UnBoundDebyeHuckelType = Potentials::UnBound::DebyeHuckel<typename Base::Topology>;
        using UnBoundDebyeHuckelSpheresType = Potentials::UnBound::DebyeHuckelSpheres<typename Base::Topology>;
        using UnBoundDebyeHuckelDistanceDependentDielectricType = Potentials::UnBound::DebyeHuckelDistanceDependentDielectric<typename Base::Topology>;
        using UnBoundDLVOType1Type = Potentials::UnBound::DLVOType1<typename Base::Topology>;
        using UnBoundDLVOType2Type = Potentials::UnBound::DLVOType2<typename Base::Topology>;
        using UnBoundDLVOType3Type = Potentials::UnBound::DLVOType3<typename Base::Topology>;
        using UnBoundClashedType = Potentials::UnBound::Clashed<typename Base::Topology>;
        using UnBoundKimHummerType = Potentials::UnBound::KimHummer<typename Base::Topology>;
        using UnBoundLennardJonesType1Type = Potentials::UnBound::LennardJonesType1<typename Base::Topology>;
        using UnBoundLennardJonesType2Type = Potentials::UnBound::LennardJonesType2<typename Base::Topology>;
        using UnBoundLennardJonesType3Type = Potentials::UnBound::LennardJonesType3<typename Base::Topology>;
        using UnBoundWCAType1Type = Potentials::UnBound::WCAType1<typename Base::Topology>;
        using UnBoundWCAType2Type = Potentials::UnBound::WCAType2<typename Base::Topology>;
        using UnBoundWCAType3Type = Potentials::UnBound::WCAType3<typename Base::Topology>;
        using UnBoundGeneralLennardJonesType1Type = Potentials::UnBound::GeneralLennardJonesType1<typename Base::Topology>;
        using UnBoundGeneralLennardJonesType2Type = Potentials::UnBound::GeneralLennardJonesType2<typename Base::Topology>;
        using UnBoundGeneralLennardJonesType3Type = Potentials::UnBound::GeneralLennardJonesType3<typename Base::Topology>;
        using UnBoundSteric6Type = Potentials::UnBound::Steric6<typename Base::Topology>;
        using UnBoundSteric12Type = Potentials::UnBound::Steric12<typename Base::Topology>;


        using InteractorBond1HarmonicType   = Interactor::BondedInteractor<Bond1HarmonicType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond1HarmonicType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond1HarmonicType>>;

        using InteractorBond1HarmonicConst_r0Type   = Interactor::BondedInteractor<Bond1HarmonicConst_r0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond1HarmonicConst_r0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond1HarmonicConst_r0Type>>;

        using InteractorBond1HarmonicConst_K_r0Type   = Interactor::BondedInteractor<Bond1HarmonicConst_K_r0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond1HarmonicConst_K_r0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond1HarmonicConst_K_r0Type>>;

        using InteractorBond1HarmonicCommon_K_r0Type   = Interactor::BondedInteractor<Bond1HarmonicCommon_K_r0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond1HarmonicCommon_K_r0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond1HarmonicCommon_K_r0Type>>;

        using InteractorBond1HarmonicConstCommon_K_r0Type   = Interactor::BondedInteractor<Bond1HarmonicConstCommon_K_r0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond1HarmonicConstCommon_K_r0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond1HarmonicConstCommon_K_r0Type>>;

        using InteractorBond2FeneType   = Interactor::BondedInteractor<Bond2FeneType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2FeneType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2FeneType>>;

        using InteractorBond2FeneConst_K_R0Type   = Interactor::BondedInteractor<Bond2FeneConst_K_R0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2FeneConst_K_R0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2FeneConst_K_R0Type>>;

        using InteractorBond2FeneConst_r0_K_R0Type   = Interactor::BondedInteractor<Bond2FeneConst_r0_K_R0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2FeneConst_r0_K_R0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2FeneConst_r0_K_R0Type>>;

        using InteractorBond2DebyeHuckelType   = Interactor::BondedInteractor<Bond2DebyeHuckelType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2DebyeHuckelType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2DebyeHuckelType>>;

        using InteractorBond2HarmonicType   = Interactor::BondedInteractor<Bond2HarmonicType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2HarmonicType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2HarmonicType>>;

        using InteractorBond2HarmonicConst_KType   = Interactor::BondedInteractor<Bond2HarmonicConst_KType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2HarmonicConst_KType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2HarmonicConst_KType>>;

        using InteractorBond2HarmonicConst_K_r0Type   = Interactor::BondedInteractor<Bond2HarmonicConst_K_r0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2HarmonicConst_K_r0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2HarmonicConst_K_r0Type>>;

        using InteractorBond2Steric6Type   = Interactor::BondedInteractor<Bond2Steric6Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2Steric6Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2Steric6Type>>;

        using InteractorBond2Steric6Const_epsilon_sigmaType   = Interactor::BondedInteractor<Bond2Steric6Const_epsilon_sigmaType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2Steric6Const_epsilon_sigmaType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2Steric6Const_epsilon_sigmaType>>;

        using InteractorBond2Steric12Type   = Interactor::BondedInteractor<Bond2Steric12Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2Steric12Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2Steric12Type>>;

        using InteractorBond2Steric12Const_epsilon_sigmaType   = Interactor::BondedInteractor<Bond2Steric12Const_epsilon_sigmaType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2Steric12Const_epsilon_sigmaType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2Steric12Const_epsilon_sigmaType>>;

        using InteractorBond2MorseType   = Interactor::BondedInteractor<Bond2MorseType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2MorseType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2MorseType>>;

        using InteractorBond2MorseConst_DType   = Interactor::BondedInteractor<Bond2MorseConst_DType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2MorseConst_DType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2MorseConst_DType>>;

        using InteractorBond2MorseConst_r0_E_DType   = Interactor::BondedInteractor<Bond2MorseConst_r0_E_DType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2MorseConst_r0_E_DType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2MorseConst_r0_E_DType>>;

        using InteractorBond2MorseWCAType   = Interactor::BondedInteractor<Bond2MorseWCAType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2MorseWCAType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2MorseWCAType>>;

        using InteractorBond2LennardJonesType2Type   = Interactor::BondedInteractor<Bond2LennardJonesType2Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2LennardJonesType2Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2LennardJonesType2Type>>;

        using InteractorBond2LennardJonesType2Const_eType   = Interactor::BondedInteractor<Bond2LennardJonesType2Const_eType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2LennardJonesType2Const_eType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2LennardJonesType2Const_eType>>;

        using InteractorBond2LennardJonesType3Type   = Interactor::BondedInteractor<Bond2LennardJonesType3Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2LennardJonesType3Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2LennardJonesType3Type>>;

        using InteractorBond2LennardJonesType3Const_eType   = Interactor::BondedInteractor<Bond2LennardJonesType3Const_eType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2LennardJonesType3Const_eType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2LennardJonesType3Const_eType>>;

        using InteractorBond2LennardJonesKaranicolasBrooksType   = Interactor::BondedInteractor<Bond2LennardJonesKaranicolasBrooksType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2LennardJonesKaranicolasBrooksType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2LennardJonesKaranicolasBrooksType>>;

        using InteractorBond2LennardJonesGaussianType   = Interactor::BondedInteractor<Bond2LennardJonesGaussianType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2LennardJonesGaussianType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2LennardJonesGaussianType>>;

        using InteractorBond2LennardJonesGaussianConst_e_DType   = Interactor::BondedInteractor<Bond2LennardJonesGaussianConst_e_DType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2LennardJonesGaussianConst_e_DType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2LennardJonesGaussianConst_e_DType>>;

        using InteractorBond2GaussianType   = Interactor::BondedInteractor<Bond2GaussianType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2GaussianType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2GaussianType>>;

        using InteractorBond2GaussianConst_E_r0_DType   = Interactor::BondedInteractor<Bond2GaussianConst_E_r0_DType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2GaussianConst_E_r0_DType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2GaussianConst_E_r0_DType>>;

        using InteractorBond2OrientedHarmonicType   = Interactor::BondedInteractor<Bond2OrientedHarmonicType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond2OrientedHarmonicType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond2OrientedHarmonicType>>;

        using InteractorBond3BestChenHummerAngularType   = Interactor::BondedInteractor<Bond3BestChenHummerAngularType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond3BestChenHummerAngularType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond3BestChenHummerAngularType>>;

        using InteractorBond3KratkyPorodType   = Interactor::BondedInteractor<Bond3KratkyPorodType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond3KratkyPorodType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond3KratkyPorodType>>;

        using InteractorBond3KratkyPorodConst_KType   = Interactor::BondedInteractor<Bond3KratkyPorodConst_KType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond3KratkyPorodConst_KType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond3KratkyPorodConst_KType>>;

        using InteractorBond3HarmonicAngularType   = Interactor::BondedInteractor<Bond3HarmonicAngularType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond3HarmonicAngularType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond3HarmonicAngularType>>;

        using InteractorBond3HarmonicAngularConst_KType   = Interactor::BondedInteractor<Bond3HarmonicAngularConst_KType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond3HarmonicAngularConst_KType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond3HarmonicAngularConst_KType>>;

        using InteractorBond3HarmonicAngularConst_K_ang0Type   = Interactor::BondedInteractor<Bond3HarmonicAngularConst_K_ang0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond3HarmonicAngularConst_K_ang0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond3HarmonicAngularConst_K_ang0Type>>;

        using InteractorBond4Dihedral4Type   = Interactor::BondedInteractor<Bond4Dihedral4Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond4Dihedral4Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond4Dihedral4Type>>;

        using InteractorBond4DihedralType   = Interactor::BondedInteractor<Bond4DihedralType,
        Interactor::BondedInteractor_ns::BondProcessor<Bond4DihedralType>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond4DihedralType>>;

        using InteractorBond4DihedralConst_n_K_phi0Type   = Interactor::BondedInteractor<Bond4DihedralConst_n_K_phi0Type,
        Interactor::BondedInteractor_ns::BondProcessor<Bond4DihedralConst_n_K_phi0Type>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,Bond4DihedralConst_n_K_phi0Type>>;

        using InteractorUnBoundDebyeHuckelType   = Interactor::PairInteractor<UnBoundDebyeHuckelType,NeighbourList>;

        using InteractorUnBoundDebyeHuckelSpheresType   = Interactor::PairInteractor<UnBoundDebyeHuckelSpheresType,NeighbourList>;

        using InteractorUnBoundDebyeHuckelDistanceDependentDielectricType   = Interactor::PairInteractor<UnBoundDebyeHuckelDistanceDependentDielectricType,NeighbourList>;

        using InteractorUnBoundDLVOType1Type   = Interactor::PairInteractor<UnBoundDLVOType1Type,NeighbourList>;

        using InteractorUnBoundDLVOType2Type   = Interactor::PairInteractor<UnBoundDLVOType2Type,NeighbourList>;

        using InteractorUnBoundDLVOType3Type   = Interactor::PairInteractor<UnBoundDLVOType3Type,NeighbourList>;

        using InteractorUnBoundClashedType   = Interactor::PairInteractor<UnBoundClashedType,NeighbourList>;

        using InteractorUnBoundKimHummerType   = Interactor::PairInteractor<UnBoundKimHummerType,NeighbourList>;

        using InteractorUnBoundLennardJonesType1Type   = Interactor::PairInteractor<UnBoundLennardJonesType1Type,NeighbourList>;

        using InteractorUnBoundLennardJonesType2Type   = Interactor::PairInteractor<UnBoundLennardJonesType2Type,NeighbourList>;

        using InteractorUnBoundLennardJonesType3Type   = Interactor::PairInteractor<UnBoundLennardJonesType3Type,NeighbourList>;

        using InteractorUnBoundWCAType1Type   = Interactor::PairInteractor<UnBoundWCAType1Type,NeighbourList>;

        using InteractorUnBoundWCAType2Type   = Interactor::PairInteractor<UnBoundWCAType2Type,NeighbourList>;

        using InteractorUnBoundWCAType3Type   = Interactor::PairInteractor<UnBoundWCAType3Type,NeighbourList>;

        using InteractorUnBoundGeneralLennardJonesType1Type   = Interactor::PairInteractor<UnBoundGeneralLennardJonesType1Type,NeighbourList>;

        using InteractorUnBoundGeneralLennardJonesType2Type   = Interactor::PairInteractor<UnBoundGeneralLennardJonesType2Type,NeighbourList>;

        using InteractorUnBoundGeneralLennardJonesType3Type   = Interactor::PairInteractor<UnBoundGeneralLennardJonesType3Type,NeighbourList>;

        using InteractorUnBoundSteric6Type   = Interactor::PairInteractor<UnBoundSteric6Type,NeighbourList>;

        using InteractorUnBoundSteric12Type   = Interactor::PairInteractor<UnBoundSteric12Type,NeighbourList>;


        std::map<std::string,std::shared_ptr<uammd::Interactor>> interactors;

        std::map<std::string,std::shared_ptr<uammd::Interactor>> stoppedInteractors;

        bool isNeighbourListInit = false;

        void initNeighbourList(InputFile& in){

            this->sys->template log<System::MESSAGE>("[Generic] Initializing neighbour list ...");
            
            VerletListDst = std::stof(in.getOption("VerletListDst",InputFile::Required).str());

            condition = std::make_shared<Condition>(this->pd,this->top,in);
            
            typename NeighbourList::Parameters NeighbourListParam;
            
            NeighbourListParam.cutOff       = real(0.0);
            NeighbourListParam.cutOffVerlet = VerletListDst;

            nl = std::make_shared<NeighbourList>(this->pg,
                                                 condition,
                                                 NeighbourListParam);
            isNeighbourListInit = true;
        }

    public:

        Generic(std::shared_ptr<ParticleGroup> pg,
                InputFile&                     in):Base(pg,in){

                this->sys->template log<System::MESSAGE>("[Generic] Start");
                
                if(this->top->isEntryPresent("Bond1","Harmonic")){

                    auto entryInfo = this->top->getEntryInfo("Bond1","Harmonic");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Harmonic");
                        
                        typename Bond1HarmonicType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond1HarmonicType> bnd = std::make_shared<Bond1HarmonicType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond1HarmonicType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond1HarmonicType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond1","HarmonicConst_r0")){

                    auto entryInfo = this->top->getEntryInfo("Bond1","HarmonicConst_r0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicConst_r0");
                        
                        typename Bond1HarmonicConst_r0Type::Parameters bndParam;
    
                        bndParam.r0=Miscellany::str2real3(entryInfo[i].param["r0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: r0 with value %s, to interactor: Bond1::HarmonicConst_r0",                                  
                                                              entryInfo[i].param["r0"].c_str());
                    
                        std::shared_ptr<Bond1HarmonicConst_r0Type> bnd = std::make_shared<Bond1HarmonicConst_r0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond1HarmonicConst_r0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond1HarmonicConst_r0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond1","HarmonicConst_K_r0")){

                    auto entryInfo = this->top->getEntryInfo("Bond1","HarmonicConst_K_r0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicConst_K_r0");
                        
                        typename Bond1HarmonicConst_K_r0Type::Parameters bndParam;
    
                        bndParam.K=Miscellany::str2real3(entryInfo[i].param["K"],this->sys);
                        bndParam.r0=Miscellany::str2real3(entryInfo[i].param["r0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond1::HarmonicConst_K_r0",                                  
                                                              entryInfo[i].param["K"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: r0 with value %s, to interactor: Bond1::HarmonicConst_K_r0",                                  
                                                              entryInfo[i].param["r0"].c_str());
                    
                        std::shared_ptr<Bond1HarmonicConst_K_r0Type> bnd = std::make_shared<Bond1HarmonicConst_K_r0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond1HarmonicConst_K_r0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond1HarmonicConst_K_r0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond1","HarmonicCommon_K_r0")){

                    auto entryInfo = this->top->getEntryInfo("Bond1","HarmonicCommon_K_r0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicCommon_K_r0");
                        
                        typename Bond1HarmonicCommon_K_r0Type::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond1HarmonicCommon_K_r0Type> bnd = std::make_shared<Bond1HarmonicCommon_K_r0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond1HarmonicCommon_K_r0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond1HarmonicCommon_K_r0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond1","HarmonicConstCommon_K_r0")){

                    auto entryInfo = this->top->getEntryInfo("Bond1","HarmonicConstCommon_K_r0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicConstCommon_K_r0");
                        
                        typename Bond1HarmonicConstCommon_K_r0Type::Parameters bndParam;
    
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        bndParam.r0=Miscellany::str2real(entryInfo[i].param["r0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond1::HarmonicConstCommon_K_r0",                                  
                                                              entryInfo[i].param["K"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: r0 with value %s, to interactor: Bond1::HarmonicConstCommon_K_r0",                                  
                                                              entryInfo[i].param["r0"].c_str());
                    
                        std::shared_ptr<Bond1HarmonicConstCommon_K_r0Type> bnd = std::make_shared<Bond1HarmonicConstCommon_K_r0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond1HarmonicConstCommon_K_r0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond1HarmonicConstCommon_K_r0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Fene")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Fene");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Fene");
                        
                        typename Bond2FeneType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2FeneType> bnd = std::make_shared<Bond2FeneType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2FeneType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2FeneType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","FeneConst_K_R0")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","FeneConst_K_R0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : FeneConst_K_R0");
                        
                        typename Bond2FeneConst_K_R0Type::Parameters bndParam;
    
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        bndParam.R0=Miscellany::str2real(entryInfo[i].param["R0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond2::FeneConst_K_R0",                                  
                                                              entryInfo[i].param["K"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: R0 with value %s, to interactor: Bond2::FeneConst_K_R0",                                  
                                                              entryInfo[i].param["R0"].c_str());
                    
                        std::shared_ptr<Bond2FeneConst_K_R0Type> bnd = std::make_shared<Bond2FeneConst_K_R0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2FeneConst_K_R0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2FeneConst_K_R0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","FeneConst_r0_K_R0")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","FeneConst_r0_K_R0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : FeneConst_r0_K_R0");
                        
                        typename Bond2FeneConst_r0_K_R0Type::Parameters bndParam;
    
                        bndParam.r0=Miscellany::str2real(entryInfo[i].param["r0"],this->sys);
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        bndParam.R0=Miscellany::str2real(entryInfo[i].param["R0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: r0 with value %s, to interactor: Bond2::FeneConst_r0_K_R0",                                  
                                                              entryInfo[i].param["r0"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond2::FeneConst_r0_K_R0",                                  
                                                              entryInfo[i].param["K"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: R0 with value %s, to interactor: Bond2::FeneConst_r0_K_R0",                                  
                                                              entryInfo[i].param["R0"].c_str());
                    
                        std::shared_ptr<Bond2FeneConst_r0_K_R0Type> bnd = std::make_shared<Bond2FeneConst_r0_K_R0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2FeneConst_r0_K_R0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2FeneConst_r0_K_R0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","DebyeHuckel")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","DebyeHuckel");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DebyeHuckel");
                        
                        typename Bond2DebyeHuckelType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2DebyeHuckelType> bnd = std::make_shared<Bond2DebyeHuckelType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2DebyeHuckelType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2DebyeHuckelType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Harmonic")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Harmonic");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Harmonic");
                        
                        typename Bond2HarmonicType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2HarmonicType> bnd = std::make_shared<Bond2HarmonicType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2HarmonicType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2HarmonicType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","HarmonicConst_K")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","HarmonicConst_K");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicConst_K");
                        
                        typename Bond2HarmonicConst_KType::Parameters bndParam;
    
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond2::HarmonicConst_K",                                  
                                                              entryInfo[i].param["K"].c_str());
                    
                        std::shared_ptr<Bond2HarmonicConst_KType> bnd = std::make_shared<Bond2HarmonicConst_KType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2HarmonicConst_KType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2HarmonicConst_KType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","HarmonicConst_K_r0")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","HarmonicConst_K_r0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicConst_K_r0");
                        
                        typename Bond2HarmonicConst_K_r0Type::Parameters bndParam;
    
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        bndParam.r0=Miscellany::str2real(entryInfo[i].param["r0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond2::HarmonicConst_K_r0",                                  
                                                              entryInfo[i].param["K"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: r0 with value %s, to interactor: Bond2::HarmonicConst_K_r0",                                  
                                                              entryInfo[i].param["r0"].c_str());
                    
                        std::shared_ptr<Bond2HarmonicConst_K_r0Type> bnd = std::make_shared<Bond2HarmonicConst_K_r0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2HarmonicConst_K_r0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2HarmonicConst_K_r0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Steric6")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Steric6");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Steric6");
                        
                        typename Bond2Steric6Type::Parameters bndParam;
    
                        bndParam.cutOff=Miscellany::str2real(entryInfo[i].param["cutOff"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOff with value %s, to interactor: Bond2::Steric6",                                  
                                                              entryInfo[i].param["cutOff"].c_str());
                    
                        std::shared_ptr<Bond2Steric6Type> bnd = std::make_shared<Bond2Steric6Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2Steric6Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2Steric6Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Steric6Const_epsilon_sigma")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Steric6Const_epsilon_sigma");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Steric6Const_epsilon_sigma");
                        
                        typename Bond2Steric6Const_epsilon_sigmaType::Parameters bndParam;
    
                        bndParam.epsilon=Miscellany::str2real(entryInfo[i].param["epsilon"],this->sys);
                        bndParam.sigma=Miscellany::str2real(entryInfo[i].param["sigma"],this->sys);
                        bndParam.cutOff=Miscellany::str2real(entryInfo[i].param["cutOff"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: epsilon with value %s, to interactor: Bond2::Steric6Const_epsilon_sigma",                                  
                                                              entryInfo[i].param["epsilon"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: sigma with value %s, to interactor: Bond2::Steric6Const_epsilon_sigma",                                  
                                                              entryInfo[i].param["sigma"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOff with value %s, to interactor: Bond2::Steric6Const_epsilon_sigma",                                  
                                                              entryInfo[i].param["cutOff"].c_str());
                    
                        std::shared_ptr<Bond2Steric6Const_epsilon_sigmaType> bnd = std::make_shared<Bond2Steric6Const_epsilon_sigmaType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2Steric6Const_epsilon_sigmaType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2Steric6Const_epsilon_sigmaType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Steric12")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Steric12");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Steric12");
                        
                        typename Bond2Steric12Type::Parameters bndParam;
    
                        bndParam.cutOff=Miscellany::str2real(entryInfo[i].param["cutOff"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOff with value %s, to interactor: Bond2::Steric12",                                  
                                                              entryInfo[i].param["cutOff"].c_str());
                    
                        std::shared_ptr<Bond2Steric12Type> bnd = std::make_shared<Bond2Steric12Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2Steric12Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2Steric12Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Steric12Const_epsilon_sigma")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Steric12Const_epsilon_sigma");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Steric12Const_epsilon_sigma");
                        
                        typename Bond2Steric12Const_epsilon_sigmaType::Parameters bndParam;
    
                        bndParam.epsilon=Miscellany::str2real(entryInfo[i].param["epsilon"],this->sys);
                        bndParam.sigma=Miscellany::str2real(entryInfo[i].param["sigma"],this->sys);
                        bndParam.cutOff=Miscellany::str2real(entryInfo[i].param["cutOff"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: epsilon with value %s, to interactor: Bond2::Steric12Const_epsilon_sigma",                                  
                                                              entryInfo[i].param["epsilon"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: sigma with value %s, to interactor: Bond2::Steric12Const_epsilon_sigma",                                  
                                                              entryInfo[i].param["sigma"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOff with value %s, to interactor: Bond2::Steric12Const_epsilon_sigma",                                  
                                                              entryInfo[i].param["cutOff"].c_str());
                    
                        std::shared_ptr<Bond2Steric12Const_epsilon_sigmaType> bnd = std::make_shared<Bond2Steric12Const_epsilon_sigmaType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2Steric12Const_epsilon_sigmaType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2Steric12Const_epsilon_sigmaType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Morse")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Morse");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Morse");
                        
                        typename Bond2MorseType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2MorseType> bnd = std::make_shared<Bond2MorseType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2MorseType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2MorseType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","MorseConst_D")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","MorseConst_D");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : MorseConst_D");
                        
                        typename Bond2MorseConst_DType::Parameters bndParam;
    
                        bndParam.D=Miscellany::str2real(entryInfo[i].param["D"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: D with value %s, to interactor: Bond2::MorseConst_D",                                  
                                                              entryInfo[i].param["D"].c_str());
                    
                        std::shared_ptr<Bond2MorseConst_DType> bnd = std::make_shared<Bond2MorseConst_DType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2MorseConst_DType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2MorseConst_DType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","MorseConst_r0_E_D")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","MorseConst_r0_E_D");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : MorseConst_r0_E_D");
                        
                        typename Bond2MorseConst_r0_E_DType::Parameters bndParam;
    
                        bndParam.r0=Miscellany::str2real(entryInfo[i].param["r0"],this->sys);
                        bndParam.E=Miscellany::str2real(entryInfo[i].param["E"],this->sys);
                        bndParam.D=Miscellany::str2real(entryInfo[i].param["D"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: r0 with value %s, to interactor: Bond2::MorseConst_r0_E_D",                                  
                                                              entryInfo[i].param["r0"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: E with value %s, to interactor: Bond2::MorseConst_r0_E_D",                                  
                                                              entryInfo[i].param["E"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: D with value %s, to interactor: Bond2::MorseConst_r0_E_D",                                  
                                                              entryInfo[i].param["D"].c_str());
                    
                        std::shared_ptr<Bond2MorseConst_r0_E_DType> bnd = std::make_shared<Bond2MorseConst_r0_E_DType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2MorseConst_r0_E_DType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2MorseConst_r0_E_DType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","MorseWCA")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","MorseWCA");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : MorseWCA");
                        
                        typename Bond2MorseWCAType::Parameters bndParam;
    
                        bndParam.eps0=Miscellany::str2real(entryInfo[i].param["eps0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: eps0 with value %s, to interactor: Bond2::MorseWCA",                                  
                                                              entryInfo[i].param["eps0"].c_str());
                    
                        std::shared_ptr<Bond2MorseWCAType> bnd = std::make_shared<Bond2MorseWCAType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2MorseWCAType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2MorseWCAType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","LennardJonesType2")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","LennardJonesType2");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesType2");
                        
                        typename Bond2LennardJonesType2Type::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2LennardJonesType2Type> bnd = std::make_shared<Bond2LennardJonesType2Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2LennardJonesType2Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2LennardJonesType2Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","LennardJonesType2Const_e")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","LennardJonesType2Const_e");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesType2Const_e");
                        
                        typename Bond2LennardJonesType2Const_eType::Parameters bndParam;
    
                        bndParam.epsilon=Miscellany::str2real(entryInfo[i].param["epsilon"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: epsilon with value %s, to interactor: Bond2::LennardJonesType2Const_e",                                  
                                                              entryInfo[i].param["epsilon"].c_str());
                    
                        std::shared_ptr<Bond2LennardJonesType2Const_eType> bnd = std::make_shared<Bond2LennardJonesType2Const_eType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2LennardJonesType2Const_eType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2LennardJonesType2Const_eType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","LennardJonesType3")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","LennardJonesType3");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesType3");
                        
                        typename Bond2LennardJonesType3Type::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2LennardJonesType3Type> bnd = std::make_shared<Bond2LennardJonesType3Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2LennardJonesType3Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2LennardJonesType3Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","LennardJonesType3Const_e")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","LennardJonesType3Const_e");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesType3Const_e");
                        
                        typename Bond2LennardJonesType3Const_eType::Parameters bndParam;
    
                        bndParam.epsilon=Miscellany::str2real(entryInfo[i].param["epsilon"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: epsilon with value %s, to interactor: Bond2::LennardJonesType3Const_e",                                  
                                                              entryInfo[i].param["epsilon"].c_str());
                    
                        std::shared_ptr<Bond2LennardJonesType3Const_eType> bnd = std::make_shared<Bond2LennardJonesType3Const_eType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2LennardJonesType3Const_eType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2LennardJonesType3Const_eType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","LennardJonesKaranicolasBrooks")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","LennardJonesKaranicolasBrooks");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesKaranicolasBrooks");
                        
                        typename Bond2LennardJonesKaranicolasBrooksType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2LennardJonesKaranicolasBrooksType> bnd = std::make_shared<Bond2LennardJonesKaranicolasBrooksType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2LennardJonesKaranicolasBrooksType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2LennardJonesKaranicolasBrooksType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","LennardJonesGaussian")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","LennardJonesGaussian");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesGaussian");
                        
                        typename Bond2LennardJonesGaussianType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2LennardJonesGaussianType> bnd = std::make_shared<Bond2LennardJonesGaussianType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2LennardJonesGaussianType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2LennardJonesGaussianType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","LennardJonesGaussianConst_e_D")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","LennardJonesGaussianConst_e_D");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesGaussianConst_e_D");
                        
                        typename Bond2LennardJonesGaussianConst_e_DType::Parameters bndParam;
    
                        bndParam.epsilon=Miscellany::str2real(entryInfo[i].param["epsilon"],this->sys);
                        bndParam.D=Miscellany::str2real(entryInfo[i].param["D"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: epsilon with value %s, to interactor: Bond2::LennardJonesGaussianConst_e_D",                                  
                                                              entryInfo[i].param["epsilon"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: D with value %s, to interactor: Bond2::LennardJonesGaussianConst_e_D",                                  
                                                              entryInfo[i].param["D"].c_str());
                    
                        std::shared_ptr<Bond2LennardJonesGaussianConst_e_DType> bnd = std::make_shared<Bond2LennardJonesGaussianConst_e_DType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2LennardJonesGaussianConst_e_DType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2LennardJonesGaussianConst_e_DType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","Gaussian")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","Gaussian");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Gaussian");
                        
                        typename Bond2GaussianType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2GaussianType> bnd = std::make_shared<Bond2GaussianType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2GaussianType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2GaussianType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","GaussianConst_E_r0_D")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","GaussianConst_E_r0_D");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : GaussianConst_E_r0_D");
                        
                        typename Bond2GaussianConst_E_r0_DType::Parameters bndParam;
    
                        bndParam.E=Miscellany::str2real(entryInfo[i].param["E"],this->sys);
                        bndParam.r0=Miscellany::str2real(entryInfo[i].param["r0"],this->sys);
                        bndParam.D=Miscellany::str2real(entryInfo[i].param["D"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: E with value %s, to interactor: Bond2::GaussianConst_E_r0_D",                                  
                                                              entryInfo[i].param["E"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: r0 with value %s, to interactor: Bond2::GaussianConst_E_r0_D",                                  
                                                              entryInfo[i].param["r0"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: D with value %s, to interactor: Bond2::GaussianConst_E_r0_D",                                  
                                                              entryInfo[i].param["D"].c_str());
                    
                        std::shared_ptr<Bond2GaussianConst_E_r0_DType> bnd = std::make_shared<Bond2GaussianConst_E_r0_DType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2GaussianConst_E_r0_DType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2GaussianConst_E_r0_DType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond2","OrientedHarmonic")){

                    auto entryInfo = this->top->getEntryInfo("Bond2","OrientedHarmonic");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : OrientedHarmonic");
                        
                        typename Bond2OrientedHarmonicType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond2OrientedHarmonicType> bnd = std::make_shared<Bond2OrientedHarmonicType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond2OrientedHarmonicType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond2OrientedHarmonicType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond3","BestChenHummerAngular")){

                    auto entryInfo = this->top->getEntryInfo("Bond3","BestChenHummerAngular");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : BestChenHummerAngular");
                        
                        typename Bond3BestChenHummerAngularType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond3BestChenHummerAngularType> bnd = std::make_shared<Bond3BestChenHummerAngularType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond3BestChenHummerAngularType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond3BestChenHummerAngularType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond3","KratkyPorod")){

                    auto entryInfo = this->top->getEntryInfo("Bond3","KratkyPorod");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : KratkyPorod");
                        
                        typename Bond3KratkyPorodType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond3KratkyPorodType> bnd = std::make_shared<Bond3KratkyPorodType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond3KratkyPorodType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond3KratkyPorodType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond3","KratkyPorodConst_K")){

                    auto entryInfo = this->top->getEntryInfo("Bond3","KratkyPorodConst_K");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : KratkyPorodConst_K");
                        
                        typename Bond3KratkyPorodConst_KType::Parameters bndParam;
    
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond3::KratkyPorodConst_K",                                  
                                                              entryInfo[i].param["K"].c_str());
                    
                        std::shared_ptr<Bond3KratkyPorodConst_KType> bnd = std::make_shared<Bond3KratkyPorodConst_KType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond3KratkyPorodConst_KType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond3KratkyPorodConst_KType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond3","HarmonicAngular")){

                    auto entryInfo = this->top->getEntryInfo("Bond3","HarmonicAngular");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicAngular");
                        
                        typename Bond3HarmonicAngularType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond3HarmonicAngularType> bnd = std::make_shared<Bond3HarmonicAngularType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond3HarmonicAngularType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond3HarmonicAngularType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond3","HarmonicAngularConst_K")){

                    auto entryInfo = this->top->getEntryInfo("Bond3","HarmonicAngularConst_K");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicAngularConst_K");
                        
                        typename Bond3HarmonicAngularConst_KType::Parameters bndParam;
    
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond3::HarmonicAngularConst_K",                                  
                                                              entryInfo[i].param["K"].c_str());
                    
                        std::shared_ptr<Bond3HarmonicAngularConst_KType> bnd = std::make_shared<Bond3HarmonicAngularConst_KType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond3HarmonicAngularConst_KType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond3HarmonicAngularConst_KType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond3","HarmonicAngularConst_K_ang0")){

                    auto entryInfo = this->top->getEntryInfo("Bond3","HarmonicAngularConst_K_ang0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : HarmonicAngularConst_K_ang0");
                        
                        typename Bond3HarmonicAngularConst_K_ang0Type::Parameters bndParam;
    
                        bndParam.ang0=Miscellany::str2real(entryInfo[i].param["ang0"],this->sys);
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: ang0 with value %s, to interactor: Bond3::HarmonicAngularConst_K_ang0",                                  
                                                              entryInfo[i].param["ang0"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond3::HarmonicAngularConst_K_ang0",                                  
                                                              entryInfo[i].param["K"].c_str());
                    
                        std::shared_ptr<Bond3HarmonicAngularConst_K_ang0Type> bnd = std::make_shared<Bond3HarmonicAngularConst_K_ang0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond3HarmonicAngularConst_K_ang0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond3HarmonicAngularConst_K_ang0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond4","Dihedral4")){

                    auto entryInfo = this->top->getEntryInfo("Bond4","Dihedral4");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Dihedral4");
                        
                        typename Bond4Dihedral4Type::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond4Dihedral4Type> bnd = std::make_shared<Bond4Dihedral4Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond4Dihedral4Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond4Dihedral4Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond4","Dihedral")){

                    auto entryInfo = this->top->getEntryInfo("Bond4","Dihedral");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Dihedral");
                        
                        typename Bond4DihedralType::Parameters bndParam;
    
                        
                        
                        std::shared_ptr<Bond4DihedralType> bnd = std::make_shared<Bond4DihedralType>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond4DihedralType::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond4DihedralType>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("Bond4","DihedralConst_n_K_phi0")){

                    auto entryInfo = this->top->getEntryInfo("Bond4","DihedralConst_n_K_phi0");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DihedralConst_n_K_phi0");
                        
                        typename Bond4DihedralConst_n_K_phi0Type::Parameters bndParam;
    
                        bndParam.n=Miscellany::str2int(entryInfo[i].param["n"],this->sys);
                        bndParam.K=Miscellany::str2real(entryInfo[i].param["K"],this->sys);
                        bndParam.phi0=Miscellany::str2real(entryInfo[i].param["phi0"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: n with value %s, to interactor: Bond4::DihedralConst_n_K_phi0",                                  
                                                              entryInfo[i].param["n"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: K with value %s, to interactor: Bond4::DihedralConst_n_K_phi0",                                  
                                                              entryInfo[i].param["K"].c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: phi0 with value %s, to interactor: Bond4::DihedralConst_n_K_phi0",                                  
                                                              entryInfo[i].param["phi0"].c_str());
                    
                        std::shared_ptr<Bond4DihedralConst_n_K_phi0Type> bnd = std::make_shared<Bond4DihedralConst_n_K_phi0Type>(this->pd,
                                                                       bndParam);
                        
                        typename InteractorBond4DihedralConst_n_K_phi0Type::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorBond4DihedralConst_n_K_phi0Type>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","DebyeHuckel")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","DebyeHuckel");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DebyeHuckel");
                        
                        typename UnBoundDebyeHuckelType::Parameters ubndParam;
    
                        ubndParam.dielectricConstant=Miscellany::str2real(entryInfo[i].param["dielectricConstant"],this->sys);
                        ubndParam.debyeLength=Miscellany::str2real(entryInfo[i].param["debyeLength"],this->sys);
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: dielectricConstant with value %s, to interactor: UnBound::DebyeHuckel",                                  
                                                              entryInfo[i].param["dielectricConstant"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: debyeLength with value %s, to interactor: UnBound::DebyeHuckel",                                  
                                                              entryInfo[i].param["debyeLength"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::DebyeHuckel",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundDebyeHuckelType> ubnd = std::make_shared<UnBoundDebyeHuckelType>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundDebyeHuckel, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundDebyeHuckelType::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundDebyeHuckel";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundDebyeHuckelType>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","DebyeHuckelSpheres")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","DebyeHuckelSpheres");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DebyeHuckelSpheres");
                        
                        typename UnBoundDebyeHuckelSpheresType::Parameters ubndParam;
    
                        ubndParam.dielectricConstant=Miscellany::str2real(entryInfo[i].param["dielectricConstant"],this->sys);
                        ubndParam.debyeLength=Miscellany::str2real(entryInfo[i].param["debyeLength"],this->sys);
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: dielectricConstant with value %s, to interactor: UnBound::DebyeHuckelSpheres",                                  
                                                              entryInfo[i].param["dielectricConstant"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: debyeLength with value %s, to interactor: UnBound::DebyeHuckelSpheres",                                  
                                                              entryInfo[i].param["debyeLength"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::DebyeHuckelSpheres",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundDebyeHuckelSpheresType> ubnd = std::make_shared<UnBoundDebyeHuckelSpheresType>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundDebyeHuckelSpheres, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundDebyeHuckelSpheresType::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundDebyeHuckelSpheres";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundDebyeHuckelSpheresType>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","DebyeHuckelDistanceDependentDielectric")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","DebyeHuckelDistanceDependentDielectric");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DebyeHuckelDistanceDependentDielectric");
                        
                        typename UnBoundDebyeHuckelDistanceDependentDielectricType::Parameters ubndParam;
    
                        ubndParam.debyeLength=Miscellany::str2real(entryInfo[i].param["debyeLength"],this->sys);
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: debyeLength with value %s, to interactor: UnBound::DebyeHuckelDistanceDependentDielectric",                                  
                                                              entryInfo[i].param["debyeLength"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::DebyeHuckelDistanceDependentDielectric",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundDebyeHuckelDistanceDependentDielectricType> ubnd = std::make_shared<UnBoundDebyeHuckelDistanceDependentDielectricType>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundDebyeHuckelDistanceDependentDielectric, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundDebyeHuckelDistanceDependentDielectricType::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundDebyeHuckelDistanceDependentDielectric";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundDebyeHuckelDistanceDependentDielectricType>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","DLVOType1")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","DLVOType1");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DLVOType1");
                        
                        typename UnBoundDLVOType1Type::Parameters ubndParam;
    
                        ubndParam.dielectricConstant=Miscellany::str2real(entryInfo[i].param["dielectricConstant"],this->sys);
                        ubndParam.debyeLength=Miscellany::str2real(entryInfo[i].param["debyeLength"],this->sys);
                        ubndParam.cutOffDstPolar=Miscellany::str2real(entryInfo[i].param["cutOffDstPolar"],this->sys);
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDstNonPolar=Miscellany::str2real(entryInfo[i].param["cutOffDstNonPolar"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: dielectricConstant with value %s, to interactor: UnBound::DLVOType1",                                  
                                                              entryInfo[i].param["dielectricConstant"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: debyeLength with value %s, to interactor: UnBound::DLVOType1",                                  
                                                              entryInfo[i].param["debyeLength"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstPolar with value %s, to interactor: UnBound::DLVOType1",                                  
                                                              entryInfo[i].param["cutOffDstPolar"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::DLVOType1",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstNonPolar with value %s, to interactor: UnBound::DLVOType1",                                  
                                                              entryInfo[i].param["cutOffDstNonPolar"].c_str());
                        
                        std::shared_ptr<UnBoundDLVOType1Type> ubnd = std::make_shared<UnBoundDLVOType1Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundDLVOType1, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundDLVOType1Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundDLVOType1";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundDLVOType1Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","DLVOType2")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","DLVOType2");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DLVOType2");
                        
                        typename UnBoundDLVOType2Type::Parameters ubndParam;
    
                        ubndParam.dielectricConstant=Miscellany::str2real(entryInfo[i].param["dielectricConstant"],this->sys);
                        ubndParam.debyeLength=Miscellany::str2real(entryInfo[i].param["debyeLength"],this->sys);
                        ubndParam.cutOffDstPolar=Miscellany::str2real(entryInfo[i].param["cutOffDstPolar"],this->sys);
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDstNonPolar=Miscellany::str2real(entryInfo[i].param["cutOffDstNonPolar"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: dielectricConstant with value %s, to interactor: UnBound::DLVOType2",                                  
                                                              entryInfo[i].param["dielectricConstant"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: debyeLength with value %s, to interactor: UnBound::DLVOType2",                                  
                                                              entryInfo[i].param["debyeLength"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstPolar with value %s, to interactor: UnBound::DLVOType2",                                  
                                                              entryInfo[i].param["cutOffDstPolar"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::DLVOType2",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstNonPolar with value %s, to interactor: UnBound::DLVOType2",                                  
                                                              entryInfo[i].param["cutOffDstNonPolar"].c_str());
                        
                        std::shared_ptr<UnBoundDLVOType2Type> ubnd = std::make_shared<UnBoundDLVOType2Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundDLVOType2, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundDLVOType2Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundDLVOType2";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundDLVOType2Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","DLVOType3")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","DLVOType3");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : DLVOType3");
                        
                        typename UnBoundDLVOType3Type::Parameters ubndParam;
    
                        ubndParam.dielectricConstant=Miscellany::str2real(entryInfo[i].param["dielectricConstant"],this->sys);
                        ubndParam.debyeLength=Miscellany::str2real(entryInfo[i].param["debyeLength"],this->sys);
                        ubndParam.cutOffDstPolar=Miscellany::str2real(entryInfo[i].param["cutOffDstPolar"],this->sys);
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDstNonPolar=Miscellany::str2real(entryInfo[i].param["cutOffDstNonPolar"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: dielectricConstant with value %s, to interactor: UnBound::DLVOType3",                                  
                                                              entryInfo[i].param["dielectricConstant"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: debyeLength with value %s, to interactor: UnBound::DLVOType3",                                  
                                                              entryInfo[i].param["debyeLength"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstPolar with value %s, to interactor: UnBound::DLVOType3",                                  
                                                              entryInfo[i].param["cutOffDstPolar"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::DLVOType3",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstNonPolar with value %s, to interactor: UnBound::DLVOType3",                                  
                                                              entryInfo[i].param["cutOffDstNonPolar"].c_str());
                        
                        std::shared_ptr<UnBoundDLVOType3Type> ubnd = std::make_shared<UnBoundDLVOType3Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundDLVOType3, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundDLVOType3Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundDLVOType3";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundDLVOType3Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","Clashed")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","Clashed");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Clashed");
                        
                        typename UnBoundClashedType::Parameters ubndParam;
    
                        ubndParam.lambda=Miscellany::str2real(entryInfo[i].param["lambda"],this->sys);
                        ubndParam.gamma=Miscellany::str2real(entryInfo[i].param["gamma"],this->sys);
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: lambda with value %s, to interactor: UnBound::Clashed",                                  
                                                              entryInfo[i].param["lambda"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: gamma with value %s, to interactor: UnBound::Clashed",                                  
                                                              entryInfo[i].param["gamma"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::Clashed",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundClashedType> ubnd = std::make_shared<UnBoundClashedType>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundClashed, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundClashedType::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundClashed";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundClashedType>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","KimHummer")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","KimHummer");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : KimHummer");
                        
                        typename UnBoundKimHummerType::Parameters ubndParam;
    
                        ubndParam.dielectricConstant=Miscellany::str2real(entryInfo[i].param["dielectricConstant"],this->sys);
                        ubndParam.debyeLenght=Miscellany::str2real(entryInfo[i].param["debyeLenght"],this->sys);
                        ubndParam.refTemperature=Miscellany::str2real(entryInfo[i].param["refTemperature"],this->sys);
                        ubndParam.epsilon_0=Miscellany::str2real(entryInfo[i].param["epsilon_0"],this->sys);
                        ubndParam.lambda=Miscellany::str2real(entryInfo[i].param["lambda"],this->sys);
                        ubndParam.sasaModel=Miscellany::str2str(entryInfo[i].param["sasaModel"],this->sys);
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.sasaLabel=Miscellany::str2str(entryInfo[i].param["sasaLabel"],this->sys);
                        ubndParam.cutOffDstNP=Miscellany::str2real(entryInfo[i].param["cutOffDstNP"],this->sys);
                        ubndParam.cutOffDstDH=Miscellany::str2real(entryInfo[i].param["cutOffDstDH"],this->sys);
                        ubndParam.zeroEnergy=Miscellany::str2real(entryInfo[i].param["zeroEnergy"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: dielectricConstant with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["dielectricConstant"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: debyeLenght with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["debyeLenght"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: refTemperature with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["refTemperature"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: epsilon_0 with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["epsilon_0"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: lambda with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["lambda"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: sasaModel with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["sasaModel"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: sasaLabel with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["sasaLabel"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstNP with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["cutOffDstNP"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDstDH with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["cutOffDstDH"].c_str());
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: zeroEnergy with value %s, to interactor: UnBound::KimHummer",                                  
                                                              entryInfo[i].param["zeroEnergy"].c_str());
                        
                        std::shared_ptr<UnBoundKimHummerType> ubnd = std::make_shared<UnBoundKimHummerType>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundKimHummer, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundKimHummerType::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundKimHummer";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundKimHummerType>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","LennardJonesType1")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","LennardJonesType1");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesType1");
                        
                        typename UnBoundLennardJonesType1Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::LennardJonesType1",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::LennardJonesType1",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundLennardJonesType1Type> ubnd = std::make_shared<UnBoundLennardJonesType1Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundLennardJonesType1, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundLennardJonesType1Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundLennardJonesType1";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundLennardJonesType1Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","LennardJonesType2")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","LennardJonesType2");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesType2");
                        
                        typename UnBoundLennardJonesType2Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::LennardJonesType2",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::LennardJonesType2",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundLennardJonesType2Type> ubnd = std::make_shared<UnBoundLennardJonesType2Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundLennardJonesType2, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundLennardJonesType2Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundLennardJonesType2";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundLennardJonesType2Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","LennardJonesType3")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","LennardJonesType3");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : LennardJonesType3");
                        
                        typename UnBoundLennardJonesType3Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::LennardJonesType3",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::LennardJonesType3",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundLennardJonesType3Type> ubnd = std::make_shared<UnBoundLennardJonesType3Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundLennardJonesType3, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundLennardJonesType3Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundLennardJonesType3";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundLennardJonesType3Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","WCAType1")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","WCAType1");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : WCAType1");
                        
                        typename UnBoundWCAType1Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::WCAType1",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::WCAType1",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundWCAType1Type> ubnd = std::make_shared<UnBoundWCAType1Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundWCAType1, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundWCAType1Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundWCAType1";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundWCAType1Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","WCAType2")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","WCAType2");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : WCAType2");
                        
                        typename UnBoundWCAType2Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::WCAType2",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::WCAType2",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundWCAType2Type> ubnd = std::make_shared<UnBoundWCAType2Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundWCAType2, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundWCAType2Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundWCAType2";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundWCAType2Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","WCAType3")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","WCAType3");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : WCAType3");
                        
                        typename UnBoundWCAType3Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::WCAType3",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::WCAType3",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundWCAType3Type> ubnd = std::make_shared<UnBoundWCAType3Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundWCAType3, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundWCAType3Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundWCAType3";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundWCAType3Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","GeneralLennardJonesType1")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","GeneralLennardJonesType1");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : GeneralLennardJonesType1");
                        
                        typename UnBoundGeneralLennardJonesType1Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::GeneralLennardJonesType1",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::GeneralLennardJonesType1",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundGeneralLennardJonesType1Type> ubnd = std::make_shared<UnBoundGeneralLennardJonesType1Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundGeneralLennardJonesType1, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundGeneralLennardJonesType1Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundGeneralLennardJonesType1";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundGeneralLennardJonesType1Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","GeneralLennardJonesType2")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","GeneralLennardJonesType2");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : GeneralLennardJonesType2");
                        
                        typename UnBoundGeneralLennardJonesType2Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::GeneralLennardJonesType2",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::GeneralLennardJonesType2",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundGeneralLennardJonesType2Type> ubnd = std::make_shared<UnBoundGeneralLennardJonesType2Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundGeneralLennardJonesType2, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundGeneralLennardJonesType2Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundGeneralLennardJonesType2";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundGeneralLennardJonesType2Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","GeneralLennardJonesType3")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","GeneralLennardJonesType3");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : GeneralLennardJonesType3");
                        
                        typename UnBoundGeneralLennardJonesType3Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::GeneralLennardJonesType3",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::GeneralLennardJonesType3",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundGeneralLennardJonesType3Type> ubnd = std::make_shared<UnBoundGeneralLennardJonesType3Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundGeneralLennardJonesType3, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundGeneralLennardJonesType3Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundGeneralLennardJonesType3";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundGeneralLennardJonesType3Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","Steric6")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","Steric6");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Steric6");
                        
                        typename UnBoundSteric6Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::Steric6",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::Steric6",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundSteric6Type> ubnd = std::make_shared<UnBoundSteric6Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundSteric6, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundSteric6Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundSteric6";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundSteric6Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }

                if(this->top->isEntryPresent("UnBound","Steric12")){
                    
                    if(!isNeighbourListInit){
                        this->initNeighbourList(in);
                    }
                    
                    auto entryInfo = this->top->getEntryInfo("UnBound","Steric12");

                    for(uint i=0;i<entryInfo.size();i++){
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : Steric12");
                        
                        typename UnBoundSteric12Type::Parameters ubndParam;
    
                        ubndParam.label=entryInfo[i].label;
                        ubndParam.cutOffDst=Miscellany::str2real(entryInfo[i].param["cutOffDst"],this->sys);
                        
                        this->sys->template log<System::MESSAGE>("[Generic] Added parameter: label with value %s, to interactor: UnBound::Steric12",                                  
                                                              entryInfo[i].label.c_str());
                    this->sys->template log<System::MESSAGE>("[Generic] Added parameter: cutOffDst with value %s, to interactor: UnBound::Steric12",                                  
                                                              entryInfo[i].param["cutOffDst"].c_str());
                        
                        std::shared_ptr<UnBoundSteric12Type> ubnd = std::make_shared<UnBoundSteric12Type>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){
                            this->sys->template log<System::CRITICAL>("[Generic] Error in potential UnBoundSteric12, cutOffDst (%f) "
                                                                      "has to be smaller than VerletListDst (%f)",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename InteractorUnBoundSteric12Type::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = "UnBoundSteric12";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<InteractorUnBoundSteric12Type>(this->pg, 
                                                                                   interactorUbndParameters);
                    }
    
                }


                if(isNeighbourListInit){
                    this->sys->template log<System::MESSAGE>("[Generic] Neighbour list cutoff: %f, Neighbour list verlet cut off: %f",
                                                              this->nl->getCutOff(),this->nl->getCutOffVerlet());
                }
        }

            
        std::vector<std::string> getComponentsList(){
            std::vector<std::string> compList;
            
            for(auto inter : interactors){
                compList.push_back(inter.first);
            }
            return compList;
        }
            
        void sum(std::string component,Computables comp,cudaStream_t st){
            for(auto inter : interactors){
                if(inter.first == component){
                    inter.second->sum(comp,st);
                }
            }
            this->sys->template log<System::CRITICAL>("[Generic] Requested potential %s to sum. "
                                                        "But %s has not been added.",
                                                        component.c_str(),component.c_str());
        }
        
        void stop(std::string alias){
            for(auto inter : interactors){
                if(inter.first == alias){
                    stoppedInteractors[alias]=interactors[alias];
                    interactors.erase(alias);
                
                    this->sys->template log<System::MESSAGE>("[Generic] "
                                                             "Stopped interactor : %s",alias.c_str());
                    return;
                }
            }
            this->sys->template log<System::CRITICAL>("[Generic] An attempt has been made to stop the %s interactor,"
                                                        "but is not present in the list of active interactors.",
                                                        alias.c_str());
        }
        
        void resume(std::string alias){
            for(auto inter : stoppedInteractors){
                if(inter.first == alias){
                    interactors[alias]=stoppedInteractors[alias];
                    stoppedInteractors.erase(alias);

                    this->sys->template log<System::MESSAGE>("[Generic] "
                                                             "Resumed interactor : %s",alias.c_str());
                    return;
                }
            }
            this->sys->template log<System::CRITICAL>("[Generic] An attempt has been made to resume the %s interactor,"
                                                        "but is not present in the list of stopped interactors.",
                                                        alias.c_str());
        }
    
        void sum(Computables comp,cudaStream_t st) override {
            Base::sum(comp,st);
            for(auto inter : interactors){
                inter.second->sum(comp,st);
            }
        }
        
        void updateBox(Box box){
            Base::updateBox(box);
            for(auto inter : interactors){
                inter.second->updateBox(box);
            }
        }

        };

}}}}
#endif

    