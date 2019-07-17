#ifdef FIX_CLASS

FixStyle(mspin/nh, FixMspinNH)

#else 

#ifndef LMP_FIX_MSPIN_NH_H
#define LMP_FIX_MSPIN_NH_H

#include "fix_rigid_nh.h"

namespace LAMMPS_NS {

  class FixMspinNH : public FixRigidNH {
    public:
      FixMspinNH(class LAMMPS *, int, char **);
      virtual ~FixMspinNH();
      virtual int setmask();
      virtual void init();
      virtual void final_integrate();
      virtual void compute_zeeman();
      virtual void compute_dipolar();
      virtual double compute_scalar();

      double extract_zeeman_pe();
      double extract_dipolar_pe();
      double extract_distance(int, int);

    protected:
      void calculate_dipoles(int);

      double **mu;
      double **dq;
      double *qm;
      int *qmcount;

      int nsum;       // total number of rigid atoms

      int zeeman_flag;
      int dipolar_flag;
      int uniform_field;
      double dipole_cutoff;

      double alpha;   // dipole interaction scaling factor
      double beta;    // zeeman+dipolar scaling factor

      double qb2f;
      double mu_0;    // force/Ampere^2 in Real
      double bxdx, bxdy, bxdz, bydx, bydy, bydz, bzdx, bzdy, bzdz;

      double zeeman_pe, dipolar_pe;
  };
}

#endif
#endif 


