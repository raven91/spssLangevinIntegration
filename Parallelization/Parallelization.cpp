//
// Created by Nikita Kruk on 12.02.20.
//

#include "../Definitions.hpp"

#if defined(MPI_PARAMETER_SCAN) \
 || defined(MPI_FAST_INTERACTION) || defined(MPI_FAST_INTERACTION_SHARED_MEMRORY)
#include <mpi.h>
#endif

#if defined(MPI_PARAMETER_SCAN) \
 || defined(MPI_FAST_INTERACTION) || defined(MPI_FAST_INTERACTION_SHARED_MEMRORY)
void LaunchParallelSession(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
}
#else
void LaunchParallelSession(int argc, char **argv)
{

}
#endif

#if defined(MPI_PARAMETER_SCAN) \
 || defined(MPI_FAST_INTERACTION) || defined(MPI_FAST_INTERACTION_SHARED_MEMRORY)
void FinalizeParallelSession()
{
  MPI_Finalize();
}
#else
void FinalizeParallelSession()
{

}
#endif