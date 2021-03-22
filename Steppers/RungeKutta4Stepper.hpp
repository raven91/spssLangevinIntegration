//
// Created by Nikita Kruk on 06.03.18.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ParticleSystemFirstOrder.hpp"
#include "../DynamicalSystems/ParticleSystemSecondOrder.hpp"
#include "../Parallelization/Thread.hpp"

#include <vector>

class RungeKutta4Stepper
{
 public:

  explicit RungeKutta4Stepper(Thread *thread, int number_of_particles, int number_of_state_variables);
  ~RungeKutta4Stepper();

  void DoStep(ParticleSystemFirstOrder &particle_system,
              std::vector<Real> &system_state,
              Real t,
              Real dt);
  void DoStep(ParticleSystemSecondOrder &particle_system,
              std::vector<Real> &system_state,
              Real t,
              Real dt);

 private:

  Thread *thread_;
  std::vector<Real> k_1_;
  std::vector<Real> k_2_;
  std::vector<Real> k_3_;
  std::vector<Real> k_4_;
  int n_;
  int s_;
  std::vector<Real> total_displacement_;

};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP
