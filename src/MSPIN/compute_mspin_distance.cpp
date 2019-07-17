#include <cstdio>
#include <cstring>
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "fix_mspin_nh.h"
#include "compute_mspin_distance.h"


using namespace LAMMPS_NS;
using namespace std;

ComputeMSDistance::ComputeMSDistance(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), rfix(NULL)
{
  if (narg != 6) error->all(FLERR, "Illegal compute rigid/mspin/distance command.");

  vector_flag = 1;
  size_vector = 1;
  extvector = 0;

  // get the fix id
  int n = strlen(arg[3]) + 1;
  rfix = new char[n];
  strcpy(rfix, arg[3]);

  ibody = force->inumeric(FLERR, arg[4]) - 1;
  jbody = force->inumeric(FLERR, arg[5]) - 1;

  memory->create(vector, size_vector, "compute/mspin:distance");
}

ComputeMSDistance::~ComputeMSDistance()
{
  delete [] rfix;
  memory->destroy(vector);
}

void ComputeMSDistance::init()
{
  irfix = modify->find_fix(rfix);
  if (irfix < 0)
    error->all(FLERR,"Fix ID for compute mspin does not exist");

  if (strncmp(modify->fix[irfix]->style,"rigid/mspin", 11))
    error->all(FLERR,"Compute mspin with non-mspin fix-ID");
}

void ComputeMSDistance::compute_vector()
{
  double dist;
  invoked_vector = update->ntimestep;

  if (strncmp(modify->fix[irfix]->style,"rigid/mspin", 11) == 0) {
    dist = ((FixMspinNH *) modify->fix[irfix])->extract_distance(ibody, jbody);
  }

  vector[0] = dist;
}