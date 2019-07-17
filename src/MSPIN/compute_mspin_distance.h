#ifdef COMPUTE_CLASS

ComputeStyle(rigid/mspin/distance, ComputeMSDistance)

#else 

#ifndef LMP_COMPUTE_MSPIN_DISTANCE_H
#define LMP_COMPUTE_MSPIN_DISTANCE_H

#include "compute.h"

namespace LAMMPS_NS {
    class ComputeMSDistance : public Compute {
        public:
            ComputeMSDistance(class LAMMPS *, int, char **);
            ~ComputeMSDistance();

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