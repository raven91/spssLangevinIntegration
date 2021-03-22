//
// Created by Nikita Kruk on 16.02.20.
//

#ifndef SPSSLANGEVININTEGRATION_STOCHASTICEULERSTEPPER_HPP
#define SPSSLANGEVININTEGRATION_STOCHASTICEULERSTEPPER_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ParticleSystemSecondOrder.hpp"
#include "../Parallelization/Thread.hpp"

#include <cmath>
#include <algorithm> // std::fill

class StochasticEulerStepper
{
 public:

  explicit StochasticEulerStepper(Thread *thread, int number_of_particles, int number_of_state_variables) :
      thread_(thread),
      buffer_(number_of_state_variables * number_of_particles, 0.0),
      n_(number_of_particles),
      s_(number_of_state_variables)
  {

  }

  ~StochasticEulerStepper()
  {
    buffer_.clear();
  }

  void DoStep(ParticleSystemSecondOrder &particle_system, std::vector<Real> &x, Real t, Real dt)
  {
    static std::vector<Real> deterministic(x.size(), 0.0), stochastic(x.size(), 0.0);
    std::fill(deterministic.begin(), deterministic.end(), 0.0);
    std::fill(stochastic.begin(), stochastic.end(), 0.0);

    particle_system.EvaluateRhs(x, deterministic, deterministic, 0.0, dt);
    particle_system.AddNoise(x, stochastic, t);

    const std::vector<int> &loop_indices = thread_->GetLoopIndices();
    for (int i : loop_indices)
    {
      int ii = s_ * i;
      for (int s = 0; s < s_; ++s)
      {
        x[ii + s] += dt * deterministic[ii + s] + std::sqrt(dt) * stochastic[ii + s];
      } // s
    } // i
  }

 private:

  Thread *thread_;
  std::vector<Real> buffer_;
  int n_;
  int s_;

};

#endif //SPSSLANGEVININTEGRATION_STOCHASTICEULERSTEPPER_HPP
