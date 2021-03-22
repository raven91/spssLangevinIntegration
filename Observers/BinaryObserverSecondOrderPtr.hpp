//
// Created by Nikita Kruk on 14.02.20.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_BINARYOBSERVERPTR_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_BINARYOBSERVERPTR_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

#include <string>
#include <fstream>
#include <chrono>

class BinaryObserverSecondOrderPtr
{
 public:

  explicit BinaryObserverSecondOrderPtr(ThreadSharedMemory *thread,
                                        Real v_0,
                                        Real xi_r,
                                        Real D_r,
                                        Real xi_phi,
                                        Real D_phi,
                                        Real sigma,
                                        Real rho,
                                        Real alpha,
                                        PeriodicBoundaryConditions &pbc_config,
                                        int number_of_particles,
                                        int number_of_state_variables,
                                        Real dt,
                                        int trial);
  ~BinaryObserverSecondOrderPtr();

  void SaveSystemState(Real *const system_state, long size, Real t);
  void SaveSummaryStatistics(const Real *const system_state, long size, Real t);

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions &pbc_config_;
  int n_;
  int s_;
  std::chrono::time_point<std::chrono::system_clock> integration_step_timer_;
  int output_time_counter_[2];
  int output_time_threshold_[2];

  std::string simulation_file_name_;
  std::ofstream simulation_file_;

  std::string summary_statistics_file_name_;
  std::ofstream summary_statistics_file_;

};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_BINARYOBSERVERPTR_HPP
