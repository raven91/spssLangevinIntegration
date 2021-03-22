//
// Created by Nikita Kruk on 14.02.20.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_SIMULATIONENGINEPTR_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_SIMULATIONENGINEPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"

#include <string>

class SimulationEngineSecondOrderPtr
{
 public:

  explicit SimulationEngineSecondOrderPtr(ThreadSharedMemory *thread,
                                          int number_of_particles,
                                          int number_of_state_variables);
  ~SimulationEngineSecondOrderPtr();

  void RunSimulation();

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions pbc_config_;
  Real *system_state_;
  long system_state_size_;
  int n_;
  int s_;

  void InitializeRandomSystemState();

};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_SIMULATIONENGINEPTR_HPP
