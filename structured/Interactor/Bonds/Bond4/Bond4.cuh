#ifndef __BOND4__
#define __BOND4__

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond4{

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
        const int l = id2index[bondParam.id_l];
        return BondType::energy(i,j,k,l,currentParticleIndex,computational,bondParam);
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
        const int l = id2index[bondParam.id_l];
        return make_real4(BondType::force(i,j,k,l,currentParticleIndex,computational,bondParam),0.0);
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
        const int l = id2index[bondParam.id_l];
        return BondType::forceTorque(i,j,k,l,currentParticleIndex,computational,bondParam);
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
        const int l = id2index[bondParam.id_l];

        real3 posi = make_real3(computational.pos[i]);
        real3 posj = make_real3(computational.pos[j]);
        real3 posk = make_real3(computational.pos[k]);
        real3 posl = make_real3(computational.pos[l]);

        ///////////////////////////////////////

        const real3 dij = computational.box.apply_pbc(posi - posj);
        const real3 djk = computational.box.apply_pbc(posj - posk);
        const real3 dlk = computational.box.apply_pbc(posl - posk);

        const real3 aijk = cross(dij,djk);
        const real3 ajkl = cross(dlk,djk);

        const real raijk2=dot(aijk,aijk);
        const real rajkl2=dot(ajkl,ajkl);

        const real inv_raijkl = rsqrt(raijk2*rajkl2);

        real cos_dih = dot(aijk,ajkl)*inv_raijkl;
        cos_dih=min(real( 1.0),cos_dih);
        cos_dih=max(real(-1.0),cos_dih);

        const real rjk     = sqrt(dot(djk,djk));

        real sin_dih = dot(aijk,dlk)*rjk*inv_raijkl;
        sin_dih=min(real( 1.0),sin_dih);
        sin_dih=max(real(-1.0),sin_dih);

        return BondType::energy(cos_dih,sin_dih,computational,bondParam)/real(4.0);
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
        const int l = id2index[bondParam.id_l];

        real3 posi = make_real3(computational.pos[i]);
        real3 posj = make_real3(computational.pos[j]);
        real3 posk = make_real3(computational.pos[k]);
        real3 posl = make_real3(computational.pos[l]);

        ///////////////////////////////////////

        const real3 dij = computational.box.apply_pbc(posi - posj);
        const real3 djk = computational.box.apply_pbc(posj - posk);
        const real3 dlk = computational.box.apply_pbc(posl - posk);

        const real3 aijk = cross(dij,djk);
        const real3 ajkl = cross(dlk,djk);

        const real raijk2=dot(aijk,aijk);
        const real rajkl2=dot(ajkl,ajkl);

        const real inv_raijkl = rsqrt(raijk2*rajkl2);

        real cos_dih = dot(aijk,ajkl)*inv_raijkl;
        cos_dih=min(real( 1.0),cos_dih);
        cos_dih=max(real(-1.0),cos_dih);

        const real rjk     = sqrt(dot(djk,djk));

        real sin_dih = dot(aijk,dlk)*rjk*inv_raijkl;
        sin_dih=min(real( 1.0),sin_dih);
        sin_dih=max(real(-1.0),sin_dih);

        real3 fi  = make_real3(0.0);
        real3 fjk = make_real3(0.0);
        real3 fl  = make_real3(0.0);

        //

        const real inv_raijk2=real(1.0)/raijk2;
        const real inv_rajkl2=real(1.0)/rajkl2;

        const real inv_rjk = real(1.0)/rjk;

        //

        const real dot_ijk = dot(dij,djk);
        const real dot_jkl = dot(djk,dlk);

        const real3 grad_i  =  rjk*inv_raijk2*aijk;
        const real3 grad_jk = (-dot_ijk*inv_raijk2*aijk+dot_jkl*inv_rajkl2*ajkl)*inv_rjk;
        const real3 grad_l  = -rjk*inv_rajkl2*ajkl;

        //

        real fmod = BondType::energyDerivate(cos_dih,sin_dih,computational,bondParam);

        fi  = fmod*grad_i;
        fjk = fmod*grad_jk;
        fl  = fmod*grad_l;

        ///////////////////////////////////////

        real3 f=make_real3(0.0);
        if (currentParticleIndex == i){
            f = fi;
        } else if (currentParticleIndex == j){
            f = -fi+fjk;
        } else if (currentParticleIndex == k){
            f = -fl-fjk;
        } else if (currentParticleIndex == l){
            f = fl;
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

    enum ParticleBondIndex {p_i = 0, p_j = 1, p_k = 2, p_l = 3, none = -1};

    AngularHessianTransverser_(tensor3* hessian,
			       const int* id,
			       const int* selectedId,
			       const int* id2index):hessian(hessian),
						    id(id),
						    selectedId(selectedId),
						    id2index(id2index){}

    //Makes the operation epsilon_abc*r_c where epsilon is the levi civita symbol
    inline __device__ tensor3 levi_civita_contraction(real3 r){
      tensor3 T = tensor3(real(0.0), -r.z,        r.y,
			  r.z,        real(0.0), -r.x,
			  -r.y,        r.x,        real(0.0));
      return T;
    }

    inline __device__ tensor3 computeHessianBox(real3 grad1, real3 grad2, tensor3 grad_2_outer_grad_1, real du2_d2dhi, real du_ddhi){
      return du2_d2dhi*outer(grad2, grad1)+du_ddhi*grad_2_outer_grad_1;
    }

    inline __device__ resultType zero(){return tensor3(0.0);}

    inline __device__ void accumulate(resultType& total,const resultType current){
        total+=current;
    }

    inline __device__ resultType compute(const int currentParticleIndex,
                                         const typename BondType::ComputationalData& computational,
                                         const typename BondType::BondParameters&    bondParam){

        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        const int k = id2index[bondParam.id_k];
        const int l = id2index[bondParam.id_l];

	const int selId     = selectedId[currentParticleIndex];

	ParticleBondIndex id1 = (currentParticleIndex == bondParam.id_i) ? p_i :
	                        (currentParticleIndex == bondParam.id_j) ? p_j :
	                        (currentParticleIndex == bondParam.id_k) ? p_k :
	                        (currentParticleIndex == bondParam.id_l) ? p_l : none;

	ParticleBondIndex id2 = (selId == bondParam.id_i) ? p_i :
	                        (selId == bondParam.id_j) ? p_j :
	                        (selId == bondParam.id_k) ? p_k :
	                        (selId == bondParam.id_l) ? p_l : none;

	if (id2==none) return tensor3();

        real3 posi = make_real3(computational.pos[i]);
        real3 posj = make_real3(computational.pos[j]);
        real3 posk = make_real3(computational.pos[k]);
        real3 posl = make_real3(computational.pos[l]);

        ///////////////////////////////////////

        const real3 dij = computational.box.apply_pbc(posi - posj);
        const real3 djk = computational.box.apply_pbc(posj - posk);
        const real3 dlk = computational.box.apply_pbc(posl - posk);
	const real3 dik = computational.box.apply_pbc(posi - posk);
	const real3 dlj = computational.box.apply_pbc(posl - posj);

        const real3 aijk = cross(dij,djk);
        const real3 ajkl = cross(dlk,djk);

        const real raijk2=dot(aijk,aijk);
        const real rajkl2=dot(ajkl,ajkl);

        const real inv_raijkl = rsqrt(raijk2*rajkl2);

        real cos_dih = dot(aijk,ajkl)*inv_raijkl;
        cos_dih=min(real( 1.0),cos_dih);
        cos_dih=max(real(-1.0),cos_dih);

        const real rjk2     = (dot(djk,djk));
	const real rij2     = (dot(dij,dij));
	const real rik2     = (dot(dik,dik));
	const real rlj2     = (dot(dlj,dlj));
	const real rlk2     = (dot(dlk,dlk));

	const real rjk      = sqrt(rjk2);
	const real rij      = sqrt(rij2);
	const real rik      = sqrt(rik2);

        real sin_dih = dot(aijk,dlk)*rjk*inv_raijkl;
        sin_dih=min(real( 1.0),sin_dih);
        sin_dih=max(real(-1.0),sin_dih);

        //

        const real inv_raijk2=real(1.0)/raijk2;
        const real inv_rajkl2=real(1.0)/rajkl2;

	const real inv_raijk=sqrt(inv_raijk2);
	const real inv_rajkl=sqrt(inv_rajkl2);

	const real inv_raijk3=inv_raijk2*inv_raijk;
	const real inv_rajkl3=inv_rajkl2*inv_rajkl;


        const real inv_rjk  = real(1.0)/rjk;

        //

        const real dot_ijk = dot(dij,djk);
        const real dot_jkl = dot(djk,dlk);

        const real3 grad_i  = -rjk*inv_raijk2*aijk;
        const real3 grad_jk = -(-dot_ijk*inv_raijk2*aijk+dot_jkl*inv_rajkl2*ajkl)*inv_rjk;
        const real3 grad_l  =  rjk*inv_rajkl2*ajkl;

	const real3 derivi_raijk = inv_raijk*(rjk2*dij+0.5*(rij2+rjk2-rik2)*djk);
	const real3 derivj_raijk = inv_raijk*(-rjk2*dij+rij2*djk-0.5*(rij2+rjk2-rik2)*(djk-dij));
	const real3 derivk_raijk = inv_raijk*(-rij2*djk-0.5*(rij2+rjk2-rik2)*(-djk+dik));

	const real3 derivl_rajkl = inv_rajkl*(rjk2*dlk-0.5*(rjk2+rlk2-rlj2)*(djk));
	const real3 derivj_rajkl = inv_rajkl*(rlk2*djk-0.5*(rjk2+rlk2-rlj2)*(dlk));
	const real3 derivk_rajkl = inv_rajkl*(-rjk2*dlk-rlk2*djk+0.5*(rjk2+rlk2-rlj2)*(djk+dlk));

	const real rjk_inv_raijk3_2 = real(2.0)*rjk*inv_raijk3;
	const real rjk_inv_rajkl3_2 = real(2.0)*rjk*inv_rajkl3;

	const tensor3 grad_i_outer_grad_i = rjk_inv_raijk3_2*outer(derivi_raijk, aijk)-
	                                    rjk*inv_raijk2*levi_civita_contraction(djk);

	const tensor3 grad_j_outer_grad_i = rjk_inv_raijk3_2*outer(derivj_raijk, aijk)+
	                                    rjk*inv_raijk2*levi_civita_contraction(dij+djk)-
	                                    inv_raijk2*inv_rjk*outer(djk,aijk);

	const tensor3 grad_k_outer_grad_i = rjk_inv_raijk3_2*outer(derivk_raijk, aijk)-
	                                    rjk*inv_raijk2*levi_civita_contraction(dij)+
	                                    inv_raijk2*inv_rjk*outer(djk,aijk);

	const tensor3 grad_l_outer_grad_j = -rjk_inv_rajkl3_2*outer(ajkl, derivj_rajkl)+
	                                    rjk*inv_rajkl2*levi_civita_contraction(dlk)+
	                                    inv_rajkl2*inv_rjk*outer(ajkl,djk);

	const tensor3 grad_l_outer_grad_k = -rjk_inv_rajkl3_2*outer(ajkl, derivk_rajkl)-
	                                    rjk*inv_rajkl2*levi_civita_contraction(dlj)-
	                                    inv_rajkl2*inv_rjk*outer(ajkl,djk);

	const tensor3 grad_l_outer_grad_l = rjk*inv_rajkl2*(levi_civita_contraction(djk) +
					    real(-2.0)*inv_rajkl*outer(derivl_rajkl, ajkl));

	const tensor3 grad_j_outer_grad_jk = -outer(djk,grad_jk)/rjk2-(dot(dij,djk)*(levi_civita_contraction(dik)*inv_raijk2+
					      real(2.0)*inv_raijk3*outer(derivj_raijk, aijk))-
					      inv_raijk2*outer(dij-djk,aijk)+
					      dot(djk,dlk)*(-levi_civita_contraction(dlk)*inv_rajkl2+
					      real(-2.0)*inv_rajkl3*outer(derivj_rajkl, ajkl))+
					      inv_rajkl2*outer(dlk,ajkl))*inv_rjk;

	const tensor3 grad_k_outer_grad_jk = outer(djk,grad_jk)/rjk2-(-dot(dij,djk)*(levi_civita_contraction(dij)*inv_raijk2-
					     real(2.0)*inv_raijk3*outer(derivk_raijk, aijk))+
					     inv_raijk2*outer(dij,aijk)+
					     dot(djk,dlk)*(levi_civita_contraction(dlj)*inv_rajkl2+
					     real(-2.0)*inv_rajkl3*outer(derivk_rajkl, ajkl))-
					     inv_rajkl2*outer(djk+dlk,ajkl))*inv_rjk;

        real du_ddhi   = BondType::energyDerivate(cos_dih, sin_dih, computational, bondParam);
	real du2_d2dhi = BondType::energySecondDerivate(cos_dih, sin_dih, computational, bondParam);

      	bool id2_smaller = id2 < id1;
	if (id2 < id1) {
	  ParticleBondIndex temp = id1;
	  id1 = id2;
	  id2 = temp;
	}

	tensor3 H;

	real3 grad1 = (id1 == p_i) ?  grad_i :
                      (id1 == p_j) ? -grad_i + grad_jk :
                      (id1 == p_k) ? -grad_l - grad_jk :
                                      grad_l;

	real3 grad2 = (id2 == p_i) ?  grad_i :
                      (id2 == p_j) ? -grad_i + grad_jk :
                      (id2 == p_k) ? -grad_l - grad_jk :
                                      grad_l;
	tensor3 grad_2_outer_grad_1;
	grad_2_outer_grad_1 = (id1 == p_i && id2 == p_i) ?  grad_i_outer_grad_i :
	                      (id1 == p_i && id2 == p_j) ?  grad_j_outer_grad_i :
	                      (id1 == p_i && id2 == p_k) ?  grad_k_outer_grad_i :
	                      (id1 == p_i && id2 == p_l) ?  tensor3() :
	                      (id1 == p_j && id2 == p_j) ? -grad_j_outer_grad_i+grad_j_outer_grad_jk :
	                      (id1 == p_j && id2 == p_k) ? -grad_l_outer_grad_j-grad_j_outer_grad_jk.transpose() :
	                      (id1 == p_j && id2 == p_l) ?  grad_l_outer_grad_j :
	                      (id1 == p_k && id2 == p_k) ? -grad_l_outer_grad_k.transpose()-grad_k_outer_grad_jk :
	                      (id1 == p_k && id2 == p_l) ?  grad_l_outer_grad_k :
	                                                    grad_l_outer_grad_l; //ll

	H = computeHessianBox(grad1, grad2, grad_2_outer_grad_1, du2_d2dhi, du_ddhi);
	if (id2_smaller) H = H.transpose();
        return H;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
      hessian[index_i] += quantity;
    }
  };

  template <class BondType_>
  struct PairwiseForceTransverser_{

    real4*    pairwiseForce;
    const int*  id;
    const int*  selectedId;
    const int*  id2index;

    using BondType   = BondType_;
    using resultType = real4;

    enum ParticleBondIndex {p_i = 0, p_j = 1, p_k = 2, p_l = 3, none = -1};

    PairwiseForceTransverser_(real4* pairwiseForce,
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

        const int i = id2index[bondParam.id_i];
        const int j = id2index[bondParam.id_j];
        const int k = id2index[bondParam.id_k];
        const int l = id2index[bondParam.id_l];

	const int selId     = selectedId[currentParticleIndex];

	ParticleBondIndex id1 = (currentParticleIndex == bondParam.id_i) ? p_i :
	                        (currentParticleIndex == bondParam.id_j) ? p_j :
	                        (currentParticleIndex == bondParam.id_k) ? p_k :
	                        (currentParticleIndex == bondParam.id_l) ? p_l : none;

	ParticleBondIndex id2 = (selId == bondParam.id_i) ? p_i :
	                        (selId == bondParam.id_j) ? p_j :
	                        (selId == bondParam.id_k) ? p_k :
	                        (selId == bondParam.id_l) ? p_l : none;

	if (id2==none) return real4();

        real3 posi = make_real3(computational.pos[i]);
        real3 posj = make_real3(computational.pos[j]);
        real3 posk = make_real3(computational.pos[k]);
        real3 posl = make_real3(computational.pos[l]);

        ///////////////////////////////////////

        const real3 dij = computational.box.apply_pbc(posi - posj);
        const real3 djk = computational.box.apply_pbc(posj - posk);
        const real3 dlk = computational.box.apply_pbc(posl - posk);
	const real3 dik = computational.box.apply_pbc(posi - posk);
	const real3 dil = computational.box.apply_pbc(posi - posl);
	const real3 dlj = computational.box.apply_pbc(posl - posj);

        const real3 aijk = cross(dij,djk);
        const real3 ajkl = cross(dlk,djk);

        const real raijk2=dot(aijk,aijk);
        const real rajkl2=dot(ajkl,ajkl);

        const real inv_raijkl = rsqrt(raijk2*rajkl2);

        real cos_dih = dot(aijk,ajkl)*inv_raijkl;
        cos_dih=min(real( 1.0),cos_dih);
        cos_dih=max(real(-1.0),cos_dih);

        const real rjk2     = (dot(djk,djk));
	const real rij2     = (dot(dij,dij));
	const real rik2     = (dot(dik,dik));
	const real rlj2     = (dot(dlj,dlj));
	const real rlk2     = (dot(dlk,dlk));

	const real rjk      = sqrt(rjk2);
	const real rij      = sqrt(rij2);
	const real rik      = sqrt(rik2);

        real sin_dih = dot(aijk,dlk)*rjk*inv_raijkl;
        sin_dih=min(real( 1.0),sin_dih);
        sin_dih=max(real(-1.0),sin_dih);

        //

        const real inv_raijk2=real(1.0)/raijk2;
        const real inv_rajkl2=real(1.0)/rajkl2;

	const real inv_raijk=sqrt(inv_raijk2);
	const real inv_rajkl=sqrt(inv_rajkl2);

	const real inv_raijk3=inv_raijk2*inv_raijk;
	const real inv_rajkl3=inv_rajkl2*inv_rajkl;


        const real inv_rjk  = real(1.0)/rjk;

        //

        const real dot_ijjk = dot(dij,djk);
        const real dot_jklk = dot(djk,dlk);
	const real dot_jkik = dot(djk,dik);
	const real dot_jklj = dot(djk,dlj);
	const real dot_lkij = dot(dlk,dij);
	const real dot_iklj = dot(dik,dlj);
	const real dot_ijik = dot(dij,dik);
	const real dot_lklj = dot(dlk,dlj);
	const real dot_ikjk = dot(dik,djk);

	const real3 fij = ( inv_raijkl*dot_jklk-cos_dih*dot_jkik*inv_raijk2)*dij;
	const real3 fik = (-inv_raijkl*dot_jklj+cos_dih*dot_ijjk*inv_raijk2)*dik;
	const real3 fil =  -rjk2*inv_raijkl*dil;
	const real3 fjk = -(inv_raijkl*(-dot_lkij-dot_iklj)+cos_dih*(inv_raijk2*dot_ijik+
								     inv_rajkl2*dot_lklj))*djk;
	const real3 fjl =  -( inv_raijkl*dot_ikjk - cos_dih*dot_jklk*inv_rajkl2)*dlj;
	const real3 fkl =  -(-inv_raijkl*dot_ijjk + cos_dih*dot_jklj*inv_rajkl2)*dlk;



        real du_ddhi   = BondType::energyDerivate(cos_dih, sin_dih, computational, bondParam);

      	bool id2_smaller = id2 < id1;
	if (id2 < id1) {
	  ParticleBondIndex temp = id1;
	  id1 = id2;
	  id2 = temp;
	}

	real3 grad_u = (id1 == p_i && id2 == p_j) ?  fij/sin_dih :
	               (id1 == p_i && id2 == p_k) ?  fik/sin_dih :
	               (id1 == p_i && id2 == p_l) ?  fil/sin_dih :
	               (id1 == p_j && id2 == p_k) ?  fjk/sin_dih :
	               (id1 == p_j && id2 == p_l) ?  fjl/sin_dih :
	               (id1 == p_k && id2 == p_l) ?  fkl/sin_dih :
	               real3(); //self terms

	real4 pairforce = make_real4(-grad_u*du_ddhi, real(0.0));
	if (id2_smaller) pairforce = -pairforce;
        return pairforce;
    }

    inline __device__ void set(const int& index_i,resultType& quantity){
      pairwiseForce[index_i] += quantity;
    }
  };



template<class BondType_>
class Bond4Base_ {

    public:

        ///////////////////////////

        //Number of particles in the bond type

        static constexpr int nPart = 4;

        std::vector<std::string> getParticleBondLabels(){
            std::vector<std::string> labels = {"id_i","id_j","id_k","id_l"};
            return labels;
        }

        ///////////////////////////

        struct BondType : public BondType_{
            //Bond parameters
            struct BondParameters : public BondType_::BondParameters {
                int id_i;
                int id_j;
                int id_k;
                int id_l;
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

        ComputationalData getComputationalData(const Computables& computables){
            return BondType::getComputationalData(this->gd,
                                                  this->pg,storage,computables);
        }

        template<typename T>
        BondParameters processBondParameters(std::map<std::string,T>& bondParametersMap){
            BondParameters param;

            static_cast<typename BondType_::BondParameters&>(param) = BondType_::processBondParameters(this->gd,bondParametersMap);
            param.id_i = bondParametersMap.at("id_i");
            param.id_j = bondParametersMap.at("id_j");
            param.id_k = bondParametersMap.at("id_k");
            param.id_l = bondParametersMap.at("id_l");

            return param;
        }

    protected:

        std::shared_ptr<GlobalData>    gd;
        std::shared_ptr<ParticleGroup> pg;

        std::shared_ptr<ExtendedParticleData> pd;

        StorageData storage;

    public:

        Bond4Base_(std::shared_ptr<GlobalData>    gd,
                   std::shared_ptr<ParticleGroup> pg,
                   DataEntry& data):gd(gd),
                                    pg(pg),pd(getExtendedParticleData(pg)){

            storage = BondType::getStorageData(gd,pg,data);
        }

};

template<class BondType_>
class Bond4_ : public Bond4Base_<BondType_>{

    public:

        using BondType = typename Bond4_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = EnergyTransverser_<BondType>;
        using ForceTransverser  = ForceTransverser_<BondType>;
        ///////////////////////////

        Bond4_(std::shared_ptr<GlobalData>    gd,
               std::shared_ptr<ParticleGroup> pg,
               DataEntry& data):Bond4Base_<BondType_>(gd,pg,data){}

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
class AngularBond4_ : public Bond4Base_<BondType_>{

    public:

        using BondType = typename Bond4_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = AngularEnergyTransverser_<BondType>;
        using ForceTransverser  = AngularForceTransverser_<BondType>;
        using PairwiseForceTransverser = PairwiseForceTransverser_<BondType>;
        ///////////////////////////

        AngularBond4_(std::shared_ptr<GlobalData>    gd,
                      std::shared_ptr<ParticleGroup> pg,
                      DataEntry& data):Bond4Base_<BondType_>(gd,pg,data){}

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
class Bond4Torque_ : public Bond4Base_<BondType_> {

    public:

        using BondType = typename Bond4_<BondType_>::BondType;

        ///////////////////////////

        //Transverser
        using EnergyTransverser = EnergyTransverser_<BondType>;
        using ForceTransverser  = ForceTorqueTransverser_<BondType>;

    public:

        Bond4Torque_(std::shared_ptr<GlobalData>    gd,
                     std::shared_ptr<ParticleGroup> pg,
                     DataEntry& data):Bond4Base_<BondType_>(gd,pg,data){}

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
  class Bond4Hessian_ : public AngularBond4_<BondType_> {

  public:

  using BondType = typename Bond4_<BondType_>::BondType;

    ///////////////////////////

    //Transverser
    using HessianTransverser = AngularHessianTransverser_<BondType>;

  public:

    Bond4Hessian_(std::shared_ptr<GlobalData>    gd,
		  std::shared_ptr<ParticleGroup> pg,
		  DataEntry& data):AngularBond4_<BondType_>(gd,pg,data){}

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

}}}}

#endif
