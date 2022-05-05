/* ----------------------------------------------------------------------
            Contributing author: Akhlak Mahmood (NC State)
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(qmag,AtomVecQMag);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_QMAG_H
#define LMP_ATOM_VEC_QMAG_H

#include "atom_vec.h"

namespace LAMMPS_NS {

// Similar to AtomVecFull with the additional qm field for magnetic charge
class AtomVecQMag : public AtomVec {
 public:
  AtomVecQMag(class LAMMPS *);
  ~AtomVecQMag();

  void grow_pointers();
  void pack_restart_pre(int);
  void pack_restart_post(int);
  void unpack_restart_init(int);
  void data_atom_post(int);

 private:
  int *num_bond, *num_angle, *num_dihedral, *num_improper;
  int **bond_type, **angle_type, **dihedral_type, **improper_type;
  int **nspecial;

  int any_bond_negative, any_angle_negative, any_dihedral_negative, any_improper_negative;
  int bond_per_atom, angle_per_atom, dihedral_per_atom, improper_per_atom;
  int *bond_negative, *angle_negative, *dihedral_negative, *improper_negative;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

*/
