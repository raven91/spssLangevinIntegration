//
// Created by Nikita Kruk on 14.02.20.
//

#include "RungeKutta4StepperPtr.hpp"

RungeKutta4StepperPtr::RungeKutta4StepperPtr(ThreadSharedMemory *thread,
                                             int number_of_particles,
                                             int number_of_state_variables) :
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

RungeKutta4StepperPtr::~RungeKutta4StepperPtr()
{
  k_1_.clear();
  k_2_.clear();
  k_3_.clear();
  k_4_.clear();
  total_displacement_.clear();
}

void RungeKutta4StepperPtr::DoStep(ParticleSystemFirstOrderPtr &particle_system,
                                   Real *const system_state,
                                   Real t,
                                   Real dt)
{
  particle_system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  particle_system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  particle_system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
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

void RungeKutta4StepperPtr::DoStep(ParticleSystemSecondOrderPtr &particle_system,
                                   Real *const system_state,
                                   Real t,
                                   Real dt)
{
  particle_system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  particle_system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  particle_system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
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
  // todo: calculate max displacement
}