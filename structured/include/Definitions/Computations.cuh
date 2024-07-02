namespace uammd{
namespace structured{

    //GROMACS convention
    //rij = rj - ri
    //Fij force over particle i due to j
    //Virial -real(0.5)*dot(rij,Fij)
    //Stress -real(0.5)*outer(rij,Fij)

    inline
    __device__ __host__ real computeVirial(const real3& rij,const real3& Fij){
      return real(0.5)*dot(-rij,Fij);
    }

    inline
    __device__ __host__ tensor3 computeStress(const real3& rij,const real3& Fij){
      return real(0.5)*outer(-rij,Fij);
    }

    inline
    __device__ __host__ tensor3 computeHessianRadialPotential(const real3& rij, const real& invr,
      						    const real& invr2, const real& dudr,
      						    const real& d2udr2){

      tensor3 H = tensor3(0.0);

      const real rxx = rij.x*rij.x*invr2;
      const real rxy = rij.x*rij.y*invr2;
      const real rxz = rij.x*rij.z*invr2;
      const real ryy = rij.y*rij.y*invr2;
      const real ryz = rij.y*rij.z*invr2;
      const real rzz = rij.z*rij.z*invr2;

      H.xx = invr*(real(-1.0)+rxx)*dudr-rxx*d2udr2;
      H.xy = rxy*invr*dudr - rxy*d2udr2;
      H.xz = rxz*invr*dudr - rxz*d2udr2;

      H.yx = H.xy;
      H.yy = invr*(real(-1.0)+ryy)*dudr-ryy*d2udr2;
      H.yz = ryz*invr*dudr - ryz*d2udr2;

      H.zx = H.xz;
      H.zy = H.yz;
      H.zz = invr*(real(-1.0)+rzz)*dudr-rzz*d2udr2;
      return H;
    }

    inline
    __device__ __host__ tensor3 computeGradientAngle2AngularPotential(const real3& rji, const real3& rjk,
      							    const real3& rik, const real& invrji,
      							    const real& invrjk, const real& invrik,
      							    const real& invrji2, const real& invrjk2,
      							    const real& invrik2, const real3& gradTheta1,
      							    const real3& gradTheta2,
      							    const real& sijk, const real& cijk,
      							    const int& id1, const int& id2){

      tensor3 grad2theta;

      tensor3 I = tensor3();
      I.xx = real(1.0);
      I.yy = real(1.0);
      I.zz = real(1.0);

      tensor3 gradtheta2_gradtheta1 = outer(gradTheta2,gradTheta1);
      real cotanijk = cijk/sijk;

      if (id1 == 0 and id2 == 0){ //ii
        tensor3 term1 = cotanijk*(-real(1.0)*outer(rji,rji)*invrji2*invrji2-gradtheta2_gradtheta1+I*invrji2);
        tensor3 term2 = -(outer(gradTheta1,rji)+outer(rji,gradTheta1))*invrji2;
        grad2theta = term1 + term2;

      } else if (id1 == 0 and id2 == 1){//ij, ji
        tensor3 term1 = outer(-cotanijk*gradTheta2 + invrji2*rji + invrjk2*rjk, gradTheta1);
        tensor3 term2 = I*(1-invrji/invrjk*cijk);
        tensor3 term3 = cijk*invrji*invrjk*outer(invrji2/invrjk2*rji-rjk,rji);
        tensor3 term4 = -invrji/invrjk*sijk*outer(gradTheta2,rji);

        tensor3 kk = outer(rjk, gradTheta2);

        grad2theta = term1+invrji*invrjk/sijk*(term2+term3+term4);

      } else if (id1 == 0 and id2 == 2){//ik, ki
        tensor3 term1 = cotanijk*(-real(1.0)*gradtheta2_gradtheta1+invrji2*invrjk2*outer(rjk,rji));
        tensor3 term2 = -(invrjk2*outer(rjk,gradTheta1)+invrji2*outer(gradTheta2,rji))-invrjk*invrji/sijk*I;
        grad2theta = term1 + term2;


      } else if (id1 == 1 and id2 == 1){//jj
        tensor3 term1 = -cotanijk*gradtheta2_gradtheta1;
        tensor3 term2 = (-2*(invrji*invrjk)+cijk*(invrji2+invrjk2))*I;
        tensor3 term3 = outer(rjk*invrjk2+rji*invrji2, (rji+rjk)*invrji*invrjk);
        tensor3 term4 = -2*cijk*(outer(rji,rji)*invrji2*invrji2+outer(rjk,rjk)*invrjk2*invrjk2);
        tensor3 term5 = sijk*outer(gradTheta2, rjk*invrjk2+rji*invrji2);
        grad2theta = term1 + (term2 + term3 + term4 + term5)/sijk;

      } else if (id1 == 1 and id2 == 2){//jk, kj
        tensor3 term1 = outer(gradTheta2, -cotanijk*gradTheta1 + invrji2*rji + invrjk2*rjk);
        tensor3 term2 = I*(1-invrjk/invrji*cijk);
        tensor3 term3 = cijk*invrjk*invrji*outer(rjk, invrjk2/invrji2*rjk-rji);
        tensor3 term4 = -invrjk/invrji*sijk*outer(rjk, gradTheta1);

        tensor3 kk = outer(rjk, gradTheta1);
        grad2theta = term1+invrji*invrjk/sijk*(term2+term3+term4);

      } else if (id1 == 2 and id2 == 2){//kk
        tensor3 term1 = cotanijk*(-real(1.0)*outer(rjk,rjk)*invrjk2*invrjk2-gradtheta2_gradtheta1+I*invrjk2);
        tensor3 term2 = -(outer(gradTheta1,rjk)+outer(rjk,gradTheta1))*invrjk2;
        grad2theta = term1 + term2;

      }
      return grad2theta;
    }

}}
