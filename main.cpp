#include "Definitions.hpp"
#include "ControlEngines/SimulationEngineFirstOrder.hpp"
#include "ControlEngines/SimulationEngineFirstOrderPtr.hpp"
#include "Parallelization/Thread.hpp"
#include "Parallelization/ThreadSharedMemory.hpp"
#include "Parallelization/Parallelization.hpp"

// initialization of the random number generator
std::mt19937 mersenne_twister_generator(std::random_device{}());

int main(int argc, char **argv)
{
  LaunchParallelSession(argc, argv);
  {
    //const int kS1 = 4; // (x,y,phi,omega)
    //const int kS2 = 6; // (x,y,v_x,v_y,phi,omega)
    int number_of_particles = 1000;
    int number_of_state_variables = 4;
#if defined(MPI_PARAMETER_SCAN)
    Thread thread(argc, argv, number_of_particles, number_of_state_variables);
    SimulationEngineFirstOrder engine(&thread, number_of_particles, number_of_state_variables);
    engine.RunSimulation();
#elif defined(MPI_FAST_INTERACTION)
    ThreadTwoSided thread(argc, argv);
    SimulationEngine engine(&thread);
    engine.RunSimulation();
#elif defined(MPI_FAST_INTERACTION_SHARED_MEMRORY)
    ThreadSharedMemory thread(argc, argv, number_of_particles, number_of_state_variables);
    SimulationEngineFirstOrderPtr engine(&thread, number_of_particles, number_of_state_variables);
    engine.RunSimulation();
//    engine.RunContinuationMethod();
#else
    Thread thread(argc, argv, number_of_particles, number_of_state_variables);
    SimulationEngineFirstOrder engine(&thread, number_of_particles, number_of_state_variables);
    engine.RunSimulation();
//    engine.RunContinuationMethod();
#endif
  }
  FinalizeParallelSession();

  return 0;
}