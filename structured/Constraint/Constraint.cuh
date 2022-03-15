#ifndef __CONSTRAINT__
#define __CONSTRAINT__

namespace uammd{
namespace structured{
namespace Constraint{

class Constraint{
    
    protected:

        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;

        std::string name;

    public:

    Constraint(std::shared_ptr<System>       sys,
               std::shared_ptr<ParticleData>  pd,
               std::shared_ptr<ParticleGroup> pg,
               std::string name):sys(sys),
                                 pd(pd),pg(pg),
                                 name(name){}

    ~Constraint(){}

    std::string getName(){return name;}
            
    virtual void init( cudaStream_t st) = 0;
    virtual void applyConstraint(cudaStream_t st) = 0;
};

struct ShakeBasic_{

    Box box;

    struct Parameters{
        Box box;
    };
    
    struct BondInfo{
        real r0;
    };

    ShakeBasic_(Parameters param):box(param.box){}
        
    static __host__ BondInfo readBond(std::istream &in){
        BondInfo bi;
        in>>bi.r0;
        return bi;
    }

};

struct ShakeBasicConst_r0_ : public ShakeBasic_{

    real r0;

    struct Parameters: public ShakeBasic_{
        real r0;
    };
    
    struct BondInfo{};

    ShakeBasicConst_r0_(Parameters param):ShakeBasic_(param),
                                          r0(param.r0){}
    
    static __host__ BondInfo readBond(std::istream &in){
        BondInfo bi;
        return bi;
    }

};

template<class ShakeBondType>
class ShakeBond: public ParameterUpdatable{

    public:
        
        static constexpr int nPart = 2;

        struct Parameters: public ShakeBondType::Parameters{};

        struct Bond{

            int i,j;

            typename ShakeBondType::BondInfo bondInfo;

        };
            
        static inline Bond readBond(std::shared_ptr<System> sys,
                                    std::string& line) { 
                
            std::stringstream parser(line);
        
            int i,j;
            if(!(parser>>i>>j)) {
                sys->log<System::CRITICAL>("Line unreable, %s", 
                                            line.c_str());
            }
        
            Bond bond;
            
            bond.i = i;
            bond.j = j;

            bond.bondInfo = ShakeBondType::readBond(parser);
        
            return bond;
        
        }

        static inline void registerBond(std::shared_ptr<System>        sys,
                                        std::vector<std::vector<int>>& bondsPerParticle,
                                        int bondIndex,
                                        Bond b){
            
            //i
            bondsPerParticle[b.i].push_back(bondIndex);
            //j
            bondsPerParticle[b.j].push_back(bondIndex);
        
        }

    private:
        
        std::shared_ptr<ParticleData> pd;
        
        Parameters param;

    public:

        ShakeBond(std::shared_ptr<ParticleData> pd,
                  Parameters param):pd(pd),
                                    param(param){}

        void updateBox(Box box){
            if constexpr (has_box<Parameters>::value) {
                param.box = box;
            }
        }
};

template<class BondType_      ,
         class BondProcessor_ ,
         class BondReader_    >
class SHAKE: public Constraint, 
             public ParameterUpdatableDelegate<BondType_> {


    protected:

        using BondType    = BondType_;
        using Bond        = typename BondType::Bond;
        
        using BondProcessor = BondProcessor_;
        using BondReader    = BondReader_;
            
        std::string bondName;
        std::string restraintName;
            
        using InputType     = typename BondReader::InputType;

        using BondListType  = Interactor::BondedInteractor_ns::BondList<BondType,
                                                               BondProcessor,
                                                               BondReader>;
            
        
        std::shared_ptr<BondType> bondType;
        
        std::shared_ptr<InputType> input;
            
        std::shared_ptr<BondListType> bondList;

    public:
        
        struct Parameters{
            std::string bondName;
            std::string restraintName;
        };

        SHAKE(std::shared_ptr<System>       sys,
              std::shared_ptr<ParticleData>  pd,
              std::shared_ptr<ParticleGroup> pg,
              std::shared_ptr<InputType>  input,
              std::shared_ptr<BondType> bondType,
              Parameters param                 ):Constraint(sys,pd,pg,"SHAKE"),
                                                 bondName(param.bondName),
                                                 restraintName(param.restraintName){
            this->setDelegate(bondType);
            
            bondList = std::make_shared<BondListType>(sys,pg,input,bondName);
        }

        ~SHAKE(){}

        void init(cudaStream_t st) override {}
        void applyConstraint(cudaStream_t st) override {
        }

};

}}}

#endif
