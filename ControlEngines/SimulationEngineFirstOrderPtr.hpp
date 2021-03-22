//
// Created by Nikita Kruk on 27.05.20.
//

#ifndef SPSSLANGEVININTEGRATION_SIMULATIONENGINEFIRSTORDERPTR_HPP
#define SPSSLANGEVININTEGRATION_SIMULATIONENGINEFIRSTORDERPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"

#include <string>
#include <vector>

class SimulationEngineFirstOrderPtr
{
 public:

  explicit SimulationEngineFirstOrderPtr(ThreadSharedMemory *thread,
                                         int number_of_particles,
                                         int number_of_state_variables);
  ~SimulationEngineFirstOrderPtr();

  void RunSimulation();
  void RunContinuationMethod();

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions pbc_config_;
  Real *system_state_;
  long system_state_size_;
  std::vector<Real> natural_frequencies_;
  int n_;
  int s_;

  void InitializeRandomSystemState(Real sigma, Real xi);
  void InitializeSystemStateFromPreviousSolution(const std::string &file_name);

};

#endif //SPSSLANGEVININTEGRATION_SIMULATIONENGINEFIRSTORDERPTR_HPP
