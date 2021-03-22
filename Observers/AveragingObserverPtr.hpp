//
// Created by Nikita Kruk on 01.09.20.
//

#ifndef SPSSLANGEVININTEGRATION_AVERAGINGOBSERVERPTR_HPP
#define SPSSLANGEVININTEGRATION_AVERAGINGOBSERVERPTR_HPP

#include "../Definitions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

#include <vector>

class AveragingObserverPtr
{
 public:

  AveragingObserverPtr(ThreadSharedMemory *thread, int number_of_particles, int number_of_state_variables);
  ~AveragingObserverPtr();

  void operator()(Real *const system_state, long size, Real t);
  void SaveAveragedSystemState(Real v_0,
                               Real xi,
                               Real sigma,
                               Real rho,
                               Real alpha,
                               Real D_phi,
                               int trial);

 private:

  ThreadSharedMemory *thread_;
  int n_;
  int s_;
  std::vector<Real> accumulated_system_state_;
  int number_of_accumulated_steps_;

};

#endif //SPSSLANGEVININTEGRATION_AVERAGINGOBSERVERPTR_HPP
