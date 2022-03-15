#ifndef __LINKERS_BONDED__
#define __LINKERS_BONDED__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Disordered{
    
    template<class Base_ >
    class LinkersBonded: public Base_ {
            
            using Base = Base_;

        public:
            
            using Topology  = typename Base::Topology;
            
            using Units = typename Topology::Units;
            using Types = typename Topology::Types;
            
        protected:
            
            static constexpr real K  = real(378);
            static constexpr real r0 = real(3.8);

            using BondType     = Potentials::Bond2::HarmonicConst_K_r0;
            using AngleType    = Potentials::Bond3::BestChenHummerAngular;
            using DihedralType = Potentials::Bond4::Dihedral4;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;
            using InteractorDihedralType   = Interactor::BondedInteractor<DihedralType,
                                                                          Interactor::BondedInteractor_ns::BondProcessor<DihedralType>,
                                                                          Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,DihedralType>>;

            std::shared_ptr<InteractorBondType>     bonds;
            std::shared_ptr<InteractorAngleType>    angles;
            std::shared_ptr<InteractorDihedralType> dihedrals;
            
        public:

            LinkersBonded(std::shared_ptr<System>       sys, 
                          std::shared_ptr<ParticleData>  pd,
                          std::shared_ptr<ParticleGroup> pg,
                          InputFile& in):Base(sys,pd,pg,in){

                //Linkers
                BondType::Parameters bondParameters;
                
                bondParameters.K   = K;
                bondParameters.r0  = r0;
                
                std::shared_ptr<BondType> bH_PBC = std::make_shared<BondType>(this->pd,
                                                                              bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "LINKER_BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd,
                                                             this->top, bH_PBC,
                                                             interactorBondParameters);
                
                //
                
                AngleType::Parameters angleParameters;
                
                std::shared_ptr<AngleType> bA_PBC = std::make_shared<AngleType>(this->pd,
                                                                                angleParameters);
                
                typename InteractorAngleType::Parameters interactorAngleParameters;
                
                interactorAngleParameters.bondName = "LINKER_ANGLES";

                angles = std::make_shared<InteractorAngleType>(this->sys, this->pd,
                                                               this->top, bA_PBC,
                                                               interactorAngleParameters);
                
                //
                
                DihedralType::Parameters dihedralParameters;
                
                std::shared_ptr<DihedralType> bD_PBC = std::make_shared<DihedralType>(this->pd,
                                                                                      dihedralParameters);
                
                typename InteractorDihedralType::Parameters interactorDihedralParameters;
                
                interactorDihedralParameters.bondName = "LINKER_DIHEDRALS";

                dihedrals = std::make_shared<InteractorDihedralType>(this->sys, this->pd,
                                                                     this->top, bD_PBC,
                                                                     interactorDihedralParameters);
                
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                bonds->sum(comp,st);
                angles->sum(comp,st);
                dihedrals->sum(comp,st);
                Base::sum(comp,st);
            }
            
            void updateBox(Box box){
                bonds->updateBox(box);
                angles->updateBox(box);
                dihedrals->updateBox(box);
                Base::updateBox(box);
            }
    };

}}}}


#endif
