//
// Created by Nikita Kruk on 14.08.20.
//

#ifndef SPSSLANGEVININTEGRATION_STOCHASTICRUNGEKUTTASTEPPERPTR_HPP
#define SPSSLANGEVININTEGRATION_STOCHASTICRUNGEKUTTASTEPPERPTR_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/ParticleSystemFirstOrderPtr.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

#include <cmath>
#include <algorithm> // std::fill
#include <random>
#include <vector>

/**
 * Numerical Solution of Stochastic Differential Equations with Jumps in Finance
 * Strong Order 1.5 Taylor Scheme
 */
class StochasticRungeKuttaStepperPtr
{
 public:

  explicit StochasticRungeKuttaStepperPtr(ThreadSharedMemory *thread,
                                          int number_of_particles,
                                          int number_of_state_variables) :
      thread_(thread),
      n_(number_of_particles),
      s_(number_of_state_variables),
      total_displacement_(number_of_state_variables * number_of_particles, 0.0)
  {

  }
  ~StochasticRungeKuttaStepperPtr()
  {
    total_displacement_.clear();
  }

  void DoStep(ParticleSystemFirstOrderPtr &particle_system, Real *const system_state, Real t, Real dt)
  {
    static std::vector<Real> derivative(s_ * n_, 0.0);
    std::fill(derivative.begin(), derivative.end(), 0.0);
    particle_system.is_higher_order_sde_integration_ = true;
    particle_system.EvaluateRhs(system_state, derivative, derivative, 0.0, dt);
    const std::vector<Real> &complimentary_alignment_force = particle_system.GetComplimentaryAlignmentForce();
    Real noise = particle_system.D_phi();
    Real friction = particle_system.xi();

    std::normal_distribution<Real> norm_dist(0.0, 1.0);
    auto &loop_indices = thread_->GetLoopIndices();
    for (int i : loop_indices)
    {
      Real U_1 = norm_dist(mersenne_twister_generator), U_2 = norm_dist(mersenne_twister_generator);
      Real delta_W_i = U_1 * sqrt(dt);
      Real delta_Z_i = 0.5 * std::pow(dt, 1.5) * (U_1 + 1.0 / std::sqrt(3.0) * U_2);

      int ii = s_ * i;
      total_displacement_[ii] = derivative[ii] * dt - 0.5 * derivative[ii + 2] * derivative[ii + 1] * dt * dt;
      total_displacement_[ii + 1] = derivative[ii + 1] * dt + 0.5 * derivative[ii + 2] * derivative[ii] * dt * dt;
      total_displacement_[ii + 2] =
          derivative[ii + 2] * dt + 0.5 * derivative[ii + 3] * dt * dt + sqrt(2.0 * noise) * delta_Z_i;
      total_displacement_[ii + 3] =
          derivative[ii + 3] * dt + 0.5 * (complimentary_alignment_force[i] - friction * derivative[ii + 3]) * dt * dt
              + std::sqrt(2.0 * noise) * (delta_W_i - friction * delta_Z_i);

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

 private:

  ThreadSharedMemory *thread_;
  int n_;
  int s_;
  std::vector<Real> total_displacement_;

};

#endif //SPSSLANGEVININTEGRATION_STOCHASTICRUNGEKUTTASTEPPERPTR_HPP
