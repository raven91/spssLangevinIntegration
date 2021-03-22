//
// Created by Nikita Kruk on 16.02.20.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_STOCHASTICEULERSTEPPERPTR_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_STOCHASTICEULERSTEPPERPTR_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ParticleSystemFirstOrderPtr.hpp"
#include "../DynamicalSystems/ParticleSystemSecondOrderPtr.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

#include <cmath>
#include <algorithm> // std::fill
#include <cassert>
#include <vector>

/**
 * Numerical Solution of Stochastic Differential Equations with Jumps in Finance
 * Strong Order 0.5 Taylor Scheme
 */
class StochasticEulerStepperPtr
{
 public:

  explicit StochasticEulerStepperPtr(ThreadSharedMemory *thread, int number_of_particles, int number_of_state_variables)
      :
      thread_(thread),
      n_(number_of_particles),
      s_(number_of_state_variables),
      total_displacement_(number_of_state_variables * number_of_particles, 0.0)
  {

  }
  ~StochasticEulerStepperPtr()
  {
    total_displacement_.clear();
  }

  void DoStep(ParticleSystemFirstOrderPtr &particle_system, Real *const system_state, Real t, Real dt)
  {
    static std::vector<Real> derivative(s_ * n_, 0.0);
    std::fill(derivative.begin(), derivative.end(), 0.0);
//    static std::vector<std::vector<Real>> additional_derivative(n_, std::vector<Real>(2, 0.0));
//    std::fill(additional_derivative.begin(), additional_derivative.end(), std::vector<Real>(2, 0.0));
//    particle_system.EvaluateRhs(system_state, derivative, additional_derivative, dt);
    particle_system.EvaluateRhs(system_state, derivative, derivative, 0.0, dt);

    std::normal_distribution<Real> norm_dist(0.0, 1.0);
    Real delta_W_rx = 0.0, delta_W_ry = 0.0, delta_W_phi = 0.0;
    static const Real intensity_r = std::sqrt(2.0 * particle_system.D_phi()),
        intensity_phi = std::sqrt(2.0 * particle_system.D_phi());

    for (int i : thread_->GetLoopIndices())
    {
      delta_W_phi = std::sqrt(dt) * norm_dist(mersenne_twister_generator);

      int ii = s_ * i;
      total_displacement_[ii] = derivative[ii] * dt;
      total_displacement_[ii + 1] = derivative[ii + 1] * dt;
      total_displacement_[ii + 2] = derivative[ii + 2] * dt;
      total_displacement_[ii + 3] = derivative[ii + 3] * dt + intensity_phi * delta_W_phi;

      system_state[ii] += total_displacement_[ii];
      system_state[ii + 1] += total_displacement_[ii + 1];
      system_state[ii + 2] += total_displacement_[ii + 2];
      system_state[ii + 3] += total_displacement_[ii + 3];
    } // i
    if (particle_system.uses_verlet_list_)
    {
      particle_system.CalculateMaxDisplacement(total_displacement_);
    }
  }

  void DoStep(ParticleSystemSecondOrderPtr &particle_system, Real *const system_state, Real t, Real dt)
  {
    static std::vector<Real> derivative(s_ *n_,
    0.0);
    std::fill(derivative.begin(), derivative.end(), 0.0);
    static std::vector<std::vector<Real>> additional_derivative(n_, std::vector<Real>(2, 0.0));
    std::fill(additional_derivative.begin(), additional_derivative.end(), std::vector<Real>(2, 0.0));
//    particle_system.EvaluateRhs(system_state, derivative, additional_derivative, dt);
    particle_system.EvaluateRhs(system_state, derivative, derivative, 0.0, dt);

    std::normal_distribution<Real> norm_dist(0.0, 1.0);
    Real delta_W_rx = 0.0, delta_W_ry = 0.0, delta_W_phi = 0.0;
    static const Real intensity_r = std::sqrt(2.0 * particle_system.D_r()),
        intensity_phi = std::sqrt(2.0 * particle_system.D_phi());

    for (int i : thread_->GetLoopIndices())
    {
      delta_W_rx = std::sqrt(dt) * norm_dist(mersenne_twister_generator);
      delta_W_ry = std::sqrt(dt) * norm_dist(mersenne_twister_generator);
      delta_W_phi = std::sqrt(dt) * norm_dist(mersenne_twister_generator);

      int ii = s_ * i;
      assert(false);
      system_state[ii] += derivative[ii] * dt;
      system_state[ii + 1] += derivative[ii + 1] * dt;
      system_state[ii + 2] += derivative[ii + 2] * dt + intensity_r * delta_W_rx;
      system_state[ii + 3] += derivative[ii + 3] * dt + intensity_r * delta_W_ry;
      system_state[ii + 4] += derivative[ii + 4] * dt;
      system_state[ii + 5] += derivative[ii + 5] * dt + intensity_phi * delta_W_phi;
    } // i
    // todo: calculate max displacement
  }

 private:

  ThreadSharedMemory *thread_;
  int n_;
  int s_;
  std::vector<Real> total_displacement_;

};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_STOCHASTICEULERSTEPPERPTR_HPP
