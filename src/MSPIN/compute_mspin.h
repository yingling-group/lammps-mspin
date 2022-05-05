/* ----------------------------------------------------------------------
            Contributing author: Akhlak Mahmood (NC State)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(mspin,ComputeMspin)

#else

#ifndef LMP_COMPUTE_MSPIN_H
#define LMP_COMPUTE_MSPIN_H

#include "compute.h"

namespace LAMMPS_NS {
    class ComputeMspin : public Compute {
      public:
        ComputeMspin(class LAMMPS *, int, char **);
        ~ComputeMspin();
        void init();
        void compute_vector();

      private:
        int irfix;
        char *rfix;

        int ibody, jbody;
    };
}

#endif
#endif