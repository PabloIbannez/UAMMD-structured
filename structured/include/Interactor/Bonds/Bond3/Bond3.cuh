#pragma once

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond3{

template <class BondType_>
struct EnergyTransverser_{

    real*  energy;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = real;

    EnergyTransverser_(real*  energy,const int* id2index):energy(energy),id2index(id2index){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&   bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        const int k = id2index[bondParam.id_k];
        return BondType::energy(i,j,k,currentParticleIndex,computational,bondParam);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        energy[index_i] += quantity;
    }
};

template <class BondType_>
struct ForceTransverser_{

    real4*  force;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = real4;

    ForceTransverser_(real4*  force,const int* id2index):force(force),id2index(id2index){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        const int k = id2index[bondParam.id_k];
        return make_real4(BondType::force(i,j,k,currentParticleIndex,computational,bondParam),0.0);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i] += quantity;
    }
};

template <class BondType_>
struct ForceTorqueTransverser_{

    real4*  force;
    real4*  torque;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = ForceTorque;

    ForceTorqueTransverser_(real4*  force,real4*  torque,const int* id2index):force(force),torque(torque),id2index(id2index){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
        total.force  += current.force;
        total.torque += current.torque;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        const int k = id2index[bondParam.id_k];
        return BondType::forceTorque(i,j,k,currentParticleIndex,computational,bondParam);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i]  += quantity.force;
        torque[index_i] += quantity.torque;
    }
};

//Angular energy and force transverser

template <class BondType_>
struct AngularEnergyTransverser_{

    real*  energy;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = real;

    AngularEnergyTransverser_(real*  energy,const int* id2index):energy(energy),id2index(id2index){}

    inline __device__ resultType zero(){return real(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&   bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        const int k = id2index[bondParam.id_k];

        real3 posi = make_real3(computational.pos[i]);
        real3 posj = make_real3(computational.pos[j]);
        real3 posk = make_real3(computational.pos[k]);

        ///////////////////////////////////////

        //         i -------- j -------- k
        //             <- rji     rjk ->
        //Compute distances and vectors
        //---rji---
        const real3 rji = computational.box.apply_pbc(posi - posj);
        const real rji2 = dot(rji, rji);
        //---rkj---
        const real3 rjk = computational.box.apply_pbc(posk - posj);
        const real rjk2 = dot(rjk, rjk);

        const real inv_rjirjk = rsqrt(rji2*rjk2);

        real cijk = dot(rji, rjk)*inv_rjirjk; //cijk = cos (theta) = rji*rkj / mod(rji)*mod(rkj)
        //Cos must stay in range
        cijk = min(real( 1.0),cijk);
        cijk = max(real(-1.0),cijk);

        const real ang = acos(cijk);

        return BondType::energy(ang,computational,bondParam)/real(3.0);
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        energy[index_i] += quantity;
    }
};

template <class BondType_>
struct AngularForceTransverser_{

    real4*  force;
    const int* id2index;

    using BondType   = BondType_;
    using resultType = real4;

    AngularForceTransverser_(real4*  force,const int* id2index):force(force),id2index(id2index){}

    inline __device__ resultType zero(){return make_real4(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){
        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        const int k = id2index[bondParam.id_k];

        real3 posi = make_real3(computational.pos[i]);
        real3 posj = make_real3(computational.pos[j]);
        real3 posk = make_real3(computational.pos[k]);

        ///////////////////////////////////////

        //         i -------- j -------- k
        //             <- rji     rjk ->
        //Compute distances and vectors
        //---rji---
        const real3 rji = computational.box.apply_pbc(posi - posj);
        const real rji2 = dot(rji, rji);
        //---rkj---
        const real3 rjk = computational.box.apply_pbc(posk - posj);
        const real rjk2 = dot(rjk, rjk);

        const real inv_rjirjk = rsqrt(rji2*rjk2);

        real cijk = dot(rji, rjk)*inv_rjirjk; //cijk = cos (theta) = rji*rkj / mod(rji)*mod(rkj)
        //Cos must stay in range
        cijk = min(real( 1.0),cijk);
        cijk = max(real(-1.0),cijk);

        const real ang = acos(cijk);

        resultType result;

        real3 fi=make_real3(0.0);
        real3 fk=make_real3(0.0);

        ///////////////////////////////////////

        real fmod = BondType::energyDerivate(ang,computational,bondParam);

        const real inv_rji2 = real(1.0)/rji2;
        const real inv_rjk2 = real(1.0)/rjk2;

        real sijk = sin(ang);
        //Sin must be bigger than 0
        //sijk = max(std::numeric_limits<real>::min(),sijk);
        sijk = max(real(1e-6),sijk);
        fmod=fmod/sijk;

        const real crji = cijk*inv_rji2;
        const real crjk = cijk*inv_rjk2;

        fi = fmod*(rjk*inv_rjirjk-rji*crji);
        fk = fmod*(rji*inv_rjirjk-rjk*crjk);

        ///////////////////////////////////////

        real3 f=make_real3(0.0);
        if(currentParticleIndex==i){
            f = fi;
        }
        else if(currentParticleIndex==j){
            f = -(fi+fk);
        }
        else if(currentParticleIndex==k){
            f = fk;
        }

        return make_real4(f,real(0.0));
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
        force[index_i] += quantity;
    }
};


template <class BondType_>
struct AngularHessianTransverser_{

    tensor3*    hessian;
    const int*  id;
    const int*  selectedId;
    const int*  id2index;

    using BondType   = BondType_;
    using resultType = tensor3;

  enum ParticleBondIndex {p_i = 0, p_j = 1, p_k = 2, none = -1};
  
    AngularHessianTransverser_(tensor3* hessian,
			       const int* id,
			       const int* selectedId,
			       const int* id2index):hessian(hessian),
						    id(id),
						    selectedId(selectedId),
						    id2index(id2index){}

    inline __device__ resultType zero(){return tensor3(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
        total+=current;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){

        // xx xy xz
        // yx yy yz
        // zx zy zz

        //We compute the box ij of the hessian

      const int i = id2index[bondParam.id_i];
      const int j = id2index[bondParam.id_j];
      const int k = id2index[bondParam.id_k];

      const int selId     = selectedId[currentParticleIndex];

      tensor3 H = tensor3(0.0);

      int id1 = (currentParticleIndex == bondParam.id_i) ? p_i :
	        (currentParticleIndex == bondParam.id_j) ? p_j :
	        (currentParticleIndex == bondParam.id_k) ? p_k : none;

      int id2 = (selId == bondParam.id_i) ? p_i :
	        (selId == bondParam.id_j) ? p_j :
	        (selId == bondParam.id_k) ? p_k : none;
      
      if (id2 == none) return tensor3();
      
      real3 posi = make_real3(computational.pos[i]);
      real3 posj = make_real3(computational.pos[j]);
      real3 posk = make_real3(computational.pos[k]);
      
      ///////////////////////////////////////
      
      //         i -------- j -------- k
      //             <- rji     rjk ->
      //Compute distances and vectors
      //---rji---
      const real3 rji = computational.box.apply_pbc(posi - posj);
      const real rji2 = dot(rji, rji);
      const real invrji2 = real(1.0)/rji2;
      const real invrji  = sqrt(invrji2);
      //---rkj---
      const real3 rjk = computational.box.apply_pbc(posk - posj);
      const real rjk2 = dot(rjk, rjk);
      const real invrjk2 = real(1.0)/rjk2;
      const real invrjk  = sqrt(invrjk2);
      
      const real3 rki = computational.box.apply_pbc(posi - posk);
      const real rki2 = dot(rki, rki);
      const real invrki2 = real(1.0)/rki2;
      const real invrki  = sqrt(invrki2);
      
      const real inv_rjirjk = rsqrt(rji2*rjk2);
      
      real cijk = dot(rji, rjk)*inv_rjirjk;
      //Cos must stay in range
      cijk = min(real( 1.0),cijk);
      cijk = max(real(-1.0),cijk);
      
      const real ang = acos(cijk);
      real sijk = sin(ang);
      //Sin must be bigger than 0
      //sijk = max(std::numeric_limits<real>::min(),sijk);
      sijk = max(real(1e-6),sijk);

      const real crji = cijk*invrji2;
      const real crjk = cijk*invrjk2;

      real3 grad_i = -(rjk*inv_rjirjk-rji*crji)/sijk;
      real3 grad_k = -(rji*inv_rjirjk-rjk*crjk)/sijk;
      
      bool id2_smaller = id2 < id1;
      if (id2_smaller) {
	int temp = id1;
	id1 = id2;
	id2 = temp;
      }
     
      real3 grad1 = (id1 == p_i) ?  grad_i :
	            (id1 == p_j) ? -grad_i - grad_k :
		    grad_k;
  
      real3 grad2 = (id2 == p_i) ?  grad_i :
	            (id2 == p_j) ? -grad_i - grad_k :
                      		    grad_k;
        
      tensor3 grad_2_outer_grad_1 = computeGradientAngle2AngularPotential(rji, rjk, rki,
									  invrji, invrjk, invrki,
									  invrji2, invrjk2, invrki2,
									  grad1, grad2,
									  sijk, cijk, id1, id2);
      
      real dudtheta   = BondType::energyDerivate(ang,computational,bondParam);
      real du2dtheta2 = BondType::energySecondDerivate(ang,computational,bondParam);
      
      H = du2dtheta2*outer(grad2, grad1)+dudtheta*grad_2_outer_grad_1;
      if (id2_smaller) H = H.transpose();
      return H;
    }

  inline __device__ void set(const int& index_i,resultType& quantity){
    hessian[index_i] += quantity;
  }
};

  template <class BondType_>
  struct AngularPairwiseForceTransverser_{
    
    real4*    pairwiseForce;
    const int*  id;
    const int*  selectedId;
    const int*  id2index;
    
    using BondType   = BondType_;
    using resultType = real4;
    enum ParticleBondIndex {p_i = 0, p_j = 1, p_k = 2, none = -1};
    
    AngularPairwiseForceTransverser_(real4* pairwiseForce,
				     const int* id,
				     const int* selectedId,
				     const int* id2index):pairwiseForce(pairwiseForce),
							  id(id),
							  selectedId(selectedId),
							  id2index(id2index){}
    
    inline __device__ resultType zero(){return real4();}
    
    inline __device__ void accumulate(resultType& total,const resultType current){
      total+=current;
    }
    
    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){

        // xx xy xz
        // yx yy yz
        // zx zy zz
      
        //We compute the box ij of the hessian
      
      const int i = id2index[bondParam.id_i];
      const int j = id2index[bondParam.id_j];
      const int k = id2index[bondParam.id_k];

      const int selId     = selectedId[currentParticleIndex];

      tensor3 H = tensor3(0.0);

      int id1 = (currentParticleIndex == bondParam.id_i) ? p_i :
	        (currentParticleIndex == bondParam.id_j) ? p_j :
	        (currentParticleIndex == bondParam.id_k) ? p_k : none;

      int id2 = (selId == bondParam.id_i) ? p_i :
	        (selId == bondParam.id_j) ? p_j :
	        (selId == bondParam.id_k) ? p_k : none;
      
      if (id2 == none) return real4();
      
      real3 posi = make_real3(computational.pos[i]);
      real3 posj = make_real3(computational.pos[j]);
      real3 posk = make_real3(computational.pos[k]);
      
      ///////////////////////////////////////
      
      //         i -------- j -------- k
      //             <- rji     rjk ->
      //Compute distances and vectors
      //---rji---
      const real3 rji    = computational.box.apply_pbc(posi - posj);
      const real rji2    = dot(rji, rji);
      const real invrji  = rsqrt(rji2);
      //---rkj---
      const real3 rjk    = computational.box.apply_pbc(posk - posj);
      const real rjk2    = dot(rjk, rjk);
      const real invrjk  = rsqrt(rjk2);
      
      const real inv_rjirjk = invrjk*invrji;
      
      real cijk = dot(rji, rjk)*inv_rjirjk;
      //Cos must stay in range
      cijk = min(real( 1.0),cijk);
      cijk = max(real(-1.0),cijk);
      
      const real ang = acos(cijk);
      real sijk = sin(ang);
      //Sin must be bigger than 0
      //sijk = max(std::numeric_limits<real>::min(),sijk);
      sijk = max(real(1e-6),sijk);
      real invsijk = real(1.0)/sijk;
      
      real dtheta_drji = -(invrjk-cijk*invrji)*invrji*invsijk;
      real dtheta_drjk = -(invrji-cijk*invrjk)*invrjk*invsijk;
      real dtheta_drik =   invrji*invrjk*invsijk;
      
      bool id2_smaller = id2 < id1;
      if (id2_smaller) {
	int temp = id1;
	id1 = id2;
	id2 = temp;
      }

      real3 gradTheta = (id1 == p_i && id2 == p_j) ? dtheta_drji*rji :
	                (id1 == p_i && id2 == p_k) ? dtheta_drik*(rji-rjk) :
	                (id1 == p_j && id2 == p_k) ? -dtheta_drjk*rjk : real3();
	
      real dudtheta   = BondType::energyDerivate(ang,computational,bondParam);
            
      real4 f = make_real4(gradTheta * dudtheta * ((id2_smaller) ? -real(1.0):real(1.0)), real(0.0));
      
      return f;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
    pairwiseForce[index_i] += quantity;
    }
  };
//

template<class BondType_>
class Bond3Base_ {

    public:

        ///////////////////////////

        //Number of particles in the bond type

        static constexpr int nPart = 3;

        std::vector<std::string> getParticleBondLabels(){
            std::vector<std::string> labels = {"id_i","id_j","id_k"};
            return labels;
        }

        ///////////////////////////

        struct BondType : public BondType_{
            //Bond parameters
            struct BondParameters : public BondType_::BondParameters {
                int id_i;
                int id_j;
                int id_k;
            };
        };

        ///////////////////////////

        //Computational data
        using ComputationalData = typename BondType::ComputationalData;

        //Potential parameters
        using StorageData       = typename BondType::StorageData;

        //Bond parameters
        using BondParameters    = typename BondType::BondParameters;

        ///////////////////////////

        ComputationalData getComputationalData(const Computables& comp,
                                               const cudaStream_t& st){
            return BondType::getComputationalData(this->gd,
                                                  this->pg,storage,comp,st);
        }

        template<typename T>
        BondParameters processBondParameters(std::map<std::string,T>& bondParametersMap){
            BondParameters param;

            static_cast<typename BondType_::BondParameters&>(param) = BondType_::processBondParameters(this->gd,bondParametersMap);
            param.id_i = bondParametersMap.at("id_i");
            param.id_j = bondParametersMap.at("id_j");
            param.id_k = bondParametersMap.at("id_k");

            return param;
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        Bond3Base_(std::shared_ptr<GlobalData>    gd,
                   std::shared_ptr<ParticleGroup> pg,
                   DataEntry& data):gd(gd),
                                    pg(pg),pd(getExtendedParticleData(pg)){

            storage = BondType::getStorageData(gd,pg,data);
        }

};

template<class BondType_>
class Bond3_ : public Bond3Base_<BondType_>{

    public:

        using BondType = typename Bond3_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = EnergyTransverser_<BondType>;
        using ForceTransverser  = ForceTransverser_<BondType>;
        ///////////////////////////

        Bond3_(std::shared_ptr<GlobalData>    gd,
               std::shared_ptr<ParticleGroup> pg,
               DataEntry& data):Bond3Base_<BondType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy       = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(energy,id2index);
        }

        ForceTransverser getForceTransverser(){

            real4*  force       = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(force,id2index);
        }
};

template<class BondType_>
class AngularBond3_ : public Bond3Base_<BondType_>{

    public:

        using BondType = typename Bond3_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = AngularEnergyTransverser_<BondType>;
        using ForceTransverser  = AngularForceTransverser_<BondType>;
        using PairwiseForceTransverser = AngularPairwiseForceTransverser_<BondType>;
        ///////////////////////////

        AngularBond3_(std::shared_ptr<GlobalData>    gd,
                      std::shared_ptr<ParticleGroup> pg,
                      DataEntry& data):Bond3Base_<BondType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy       = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(energy,id2index);
        }

        ForceTransverser getForceTransverser(){

            real4*  force       = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(force,id2index);
        }

  PairwiseForceTransverser getPairwiseForceTransverser(){

            real4*  pforce        = this->pd->getPairwiseForce(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index   = this->pd->getIdOrderedIndices(access::location::gpu);
	    const int* id         = this->pd->getId(access::location::gpu, access::mode::read).raw();
            const int* selectedId = this->pd->getSelectedId(access::location::gpu, access::mode::read).raw();
            
	    
            return PairwiseForceTransverser(pforce,
					    id,
					    selectedId,id2index);
        }
};

  template<class BondType_>
  class Bond3Hessian_ : public AngularBond3_<BondType_> {

  public:

    using BondType = typename Bond3_<BondType_>::BondType;

    ///////////////////////////

    //Transverser
    using HessianTransverser = AngularHessianTransverser_<BondType>;

  public:

    Bond3Hessian_(std::shared_ptr<GlobalData>    gd,
		  std::shared_ptr<ParticleGroup> pg,
		  DataEntry& data):AngularBond3_<BondType_>(gd,pg,data){}

    ///////////////////////////

    HessianTransverser getHessianTransverser(){

      tensor3*  hessian     = this->pd->getHessian(access::location::gpu, access::mode::readwrite).raw();
      const int* id         = this->pd->getId(access::location::gpu, access::mode::read).raw();
      const int* selectedId = this->pd->getSelectedId(access::location::gpu, access::mode::read).raw();
      const int* id2index   = this->pd->getIdOrderedIndices(access::location::gpu);

      return HessianTransverser(hessian,
				id,
				selectedId,id2index);
    }
  };

  template<class BondType_>
class Bond3Torque_ : public Bond3Base_<BondType_> {

    public:

        using BondType = typename Bond3_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = EnergyTransverser_<BondType>;
        using ForceTransverser  = ForceTorqueTransverser_<BondType>;

    public:

        Bond3Torque_(std::shared_ptr<GlobalData>    gd,
                     std::shared_ptr<ParticleGroup> pg,
                     DataEntry& data):Bond3Base_<BondType_>(gd,pg,data){}

        ///////////////////////////

        EnergyTransverser getEnergyTransverser(){

            real*  energy       = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return EnergyTransverser(energy,id2index);
        }

        ForceTransverser getForceTransverser(){

            real4*  force       = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();
            const int* id2index = this->pd->getIdOrderedIndices(access::location::gpu);

            return ForceTransverser(force,id2index);
        }
};

}}}}

