//
// Created by Nikita Kruk on 14.02.20.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_RUNGEKUTTA4STEPPERPTR_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_RUNGEKUTTA4STEPPERPTR_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ParticleSystemFirstOrderPtr.hpp"
#include "../DynamicalSystems/ParticleSystemSecondOrderPtr.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

#include <vector>

class RungeKutta4StepperPtr
{
 public:

  explicit RungeKutta4StepperPtr(ThreadSharedMemory *thread, int number_of_particles, int number_of_state_variables);
  ~RungeKutta4StepperPtr();

  void DoStep(ParticleSystemFirstOrderPtr &particle_system, Real *const system_state, Real t, Real dt);
  void DoStep(ParticleSystemSecondOrderPtr &particle_system, Real *const system_state, Real t, Real dt);

 private:

  ThreadSharedMemory *thread_;
  std::vector<Real> k_1_;
  std::vector<Real> k_2_;
  std::vector<Real> k_3_;
  std::vector<Real> k_4_;
  int n_;
  int s_;
  std::vector<Real> total_displacement_;
};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_RUNGEKUTTA4STEPPERPTR_HPP
