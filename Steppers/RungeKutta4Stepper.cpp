//
// Created by Nikita Kruk on 06.03.18.
//

#include "RungeKutta4Stepper.hpp"

RungeKutta4Stepper::RungeKutta4Stepper(Thread *thread, int number_of_particles, int number_of_state_variables) :
    thread_(thread),
    k_1_(number_of_state_variables * number_of_particles, 0.0),
    k_2_(number_of_state_variables * number_of_particles, 0.0),
    k_3_(number_of_state_variables * number_of_particles, 0.0),
    k_4_(number_of_state_variables * number_of_particles, 0.0),
    n_(number_of_particles),
    s_(number_of_state_variables),
    total_displacement_(number_of_state_variables * number_of_particles, 0.0)
{

}

RungeKutta4Stepper::~RungeKutta4Stepper()
{
  k_1_.clear();
  k_2_.clear();
  k_3_.clear();
  k_4_.clear();
  total_displacement_.clear();
}

void RungeKutta4Stepper::DoStep(ParticleSystemFirstOrder &particle_system,
                                std::vector<Real> &system_state,
                                Real t,
                                Real dt)
{
  // k_1
  // No need to reset k_1_ to 0.0, since k_coef is 0.0
  particle_system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  thread_->SynchronizeVectorThoughBuffer(k_1_, k_4_);

  // k_2
  particle_system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  thread_->SynchronizeVectorThoughBuffer(k_2_, k_4_);

  // k_3
  particle_system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
  thread_->SynchronizeVectorThoughBuffer(k_3_, k_4_);

  // k_4
  particle_system.EvaluateRhs(system_state, k_3_, k_4_, 1.0, dt);

  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    int ii = s_ * i;
    for (int s = 0; s < s_; ++s)
    {
      total_displacement_[ii + s] = (k_1_[ii + s] + 2.0 * k_2_[ii + s] + 2.0 * k_3_[ii + s] + k_4_[ii + s]) * dt / 6.0;
      system_state[ii + s] += total_displacement_[ii + s];
    } // s
  } // i
  if (particle_system.uses_verlet_list_)
  {
    particle_system.CalculateMaxDisplacement(total_displacement_);
  }
}

void RungeKutta4Stepper::DoStep(ParticleSystemSecondOrder &particle_system,
                                std::vector<Real> &system_state,
                                Real t,
                                Real dt)
{
  // k_1
//	chimera_system.first_coefficient_ = true;
  // No need to reset k_1_ to 0.0, since k_coef is 0.0
//  std::fill(k_1_.begin(), k_1_.end(), 0.0);
  particle_system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  thread_->SynchronizeVectorThoughBuffer(k_1_, k_4_);
//	chimera_system.first_coefficient_ = false;

  // k_2
  particle_system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  thread_->SynchronizeVectorThoughBuffer(k_2_, k_4_);

  // k_3
  particle_system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
  thread_->SynchronizeVectorThoughBuffer(k_3_, k_4_);

  // k_4
  particle_system.last_coefficient_ = true;
  particle_system.EvaluateRhs(system_state, k_3_, k_4_, 1.0, dt);
  particle_system.last_coefficient_ = false;

  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    int ii = s_ * i;
    for (int s = 0; s < s_; ++s)
    {
      system_state[ii + s] += (k_1_[ii + s] + 2.0 * k_2_[ii + s] + 2.0 * k_3_[ii + s] + k_4_[ii + s]) * dt / 6.0;
    } // s
  } // i
}