//
// Created by Nikita Kruk on 06.03.18.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_SIMULATIONENGINE_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_SIMULATIONENGINE_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/Thread.hpp"

#include <string>

class SimulationEngineSecondOrder
{
 public:

  explicit SimulationEngineSecondOrder(Thread *thread, int number_of_particles, int number_of_state_variables);
  ~SimulationEngineSecondOrder();

  void RunSimulation();

 private:

  Thread *thread_;
  PeriodicBoundaryConditions pbc_config_;
  std::vector<Real> system_state_;
  int n_;
  int s_;

  void InitializeRandomSystemState();
  void InitializeSystemStateFromFile(const std::string &file_name);

};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_SIMULATIONENGINE_HPP
