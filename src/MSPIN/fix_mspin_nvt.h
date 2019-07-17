#ifdef FIX_CLASS

FixStyle(rigid/mspin, FixMspinNVT)

#else 

#ifndef LMP_FIX_MSPIN_NVT_H
#define LMP_FIX_MSPIN_NVT_H

#include "fix_mspin_nh.h"

namespace LAMMPS_NS
{
  class FixMspinNVT : public FixMspinNH {
    public:
      FixMspinNVT(class LAMMPS *, int, char **);
      ~FixMspinNVT() {}
  };
} // LAMMPS_NS

#endif
#endif
