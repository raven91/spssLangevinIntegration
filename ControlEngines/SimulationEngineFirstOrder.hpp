//
// Created by Nikita Kruk on 27.05.20.
//

#ifndef SPSSLANGEVININTEGRATION_SIMULATIONENGINEFIRSTORDER_HPP
#define SPSSLANGEVININTEGRATION_SIMULATIONENGINEFIRSTORDER_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/Thread.hpp"

#include <string>

class SimulationEngineFirstOrder
{
 public:

  explicit SimulationEngineFirstOrder(Thread *thread, int number_of_particles, int number_of_state_variables);
  ~SimulationEngineFirstOrder();

  void RunSimulation();
  void RunContinuationMethod();

 private:

  Thread *thread_;
  PeriodicBoundaryConditions pbc_config_;
  std::vector<Real> system_state_;
  std::vector<Real> natural_frequencies_;
  int n_;
  int s_;

  void InitializeRandomSystemState(Real sigma, Real xi);
  void InitializeSystemStateFromPreviousSolution(const std::string &file_name);

};

#endif //SPSSLANGEVININTEGRATION_SIMULATIONENGINEFIRSTORDER_HPP
