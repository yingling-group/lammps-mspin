/* ----------------------------------------------------------------------
            Contributing author: Akhlak Mahmood (NC State)
------------------------------------------------------------------------- */

#include <cstdio>
#include <cstring>
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "compute_mspin.h"
#include "fix_mspin_nh.h"

using namespace LAMMPS_NS;
using namespace std;

// @todo: change this into rigid/mspin/energy
ComputeMspin::ComputeMspin(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), rfix(NULL)
{
  if (narg != 4) error->all(FLERR, "Illegal compute mspin command.");

  vector_flag = 1;
  size_vector = 2;

  // are the quantities extensive or intensive
  extvector = 0;
  extlist = new int[size_vector];
  extlist[0] = 0;
  extlist[1] = 0;

  int n = strlen(arg[3]) + 1;
  rfix = new char[n];
  strcpy(rfix, arg[3]);

  memory->create(vector,2,"compute/mspin:vector");
}

ComputeMspin::~ComputeMspin()
{
  delete [] rfix;
  delete [] extlist;
  memory->destroy(vector);
}

void ComputeMspin::init()
{
  irfix = modify->find_fix(rfix);
  if (irfix < 0)
    error->all(FLERR,"Fix ID for compute mspin does not exist");

  if (strncmp(modify->fix[irfix]->style,"rigid/mspin",11))
    error->all(FLERR,"Compute mspin with non-mspin fix-ID");
}

void ComputeMspin::compute_vector()
{
  double zeeman, dipolar;

  invoked_vector = update->ntimestep;

  if (strncmp(modify->fix[irfix]->style,"rigid/mspin",11) == 0) {
    zeeman = ((FixMspinNH *) modify->fix[irfix])->extract_zeeman_pe();
    dipolar = ((FixMspinNH *) modify->fix[irfix])->extract_dipolar_pe();
  }

  vector[0] = dipolar;
  vector[1] = zeeman;
}