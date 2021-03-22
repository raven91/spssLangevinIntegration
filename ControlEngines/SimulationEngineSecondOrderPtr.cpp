//
// Created by Nikita Kruk on 14.02.20.
//

#include "SimulationEngineSecondOrderPtr.hpp"
#include "../Steppers/RungeKutta4StepperPtr.hpp"
#include "../Observers/BinaryObserverSecondOrderPtr.hpp"

#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm> // std::copy
#include <iostream>

SimulationEngineSecondOrderPtr::SimulationEngineSecondOrderPtr(ThreadSharedMemory *thread,
                                                               int number_of_particles,
                                                               int number_of_state_variables) :
    thread_(thread),
    pbc_config_(thread, number_of_particles, number_of_state_variables, 1.0, 1.0),
    system_state_size_(number_of_state_variables * number_of_particles),
    n_(number_of_particles),
    s_(number_of_state_variables)
{
  system_state_ = new Real[system_state_size_];
}

SimulationEngineSecondOrderPtr::~SimulationEngineSecondOrderPtr()
{
  delete[] system_state_;
}

void SimulationEngineSecondOrderPtr::InitializeRandomSystemState()
{
  if (thread_->IsRoot())
  {
    const Real two_pi = 2.0 * M_PI;

    std::uniform_real_distribution<Real> unif_real_dist_01(0.0, 1.0);
    std::uniform_real_distribution<Real> unif_real_dist_02pi(0.0, two_pi);

    for (int i = 0; i < n_; ++i)
    {
      system_state_[s_ * i] = unif_real_dist_01(mersenne_twister_generator);
      system_state_[s_ * i + 1] = unif_real_dist_01(mersenne_twister_generator);
      Real velocity_direction = unif_real_dist_02pi(mersenne_twister_generator);
      system_state_[s_ * i + 2] = std::cos(velocity_direction);
      system_state_[s_ * i + 3] = std::sin(velocity_direction);
      system_state_[s_ * i + 4] = velocity_direction;
      system_state_[s_ * i + 5] = 0.0;
    } // i

    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_, system_state_size_);
}

void SimulationEngineSecondOrderPtr::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

  Real v_0 = 1.0;
  Real xi_r = 1.0;
  Real D_r = 0.0;
  Real xi_phi = 1.0;
  Real D_phi = 0.0;
  Real sigma = 1.0;
  Real rho = 0.1;
  Real alpha = 1.5;

  Real t_0 = 0.0;
  Real t_1 = 100.0;
  Real dt = 0.01;

//  int trial = 0;
  for (int trial = 0; trial < 16; ++trial)
  {
    alpha = 0.1 * trial;
    InitializeRandomSystemState();
    RungeKutta4StepperPtr stepper(thread_, n_, s_);
    ParticleSystemSecondOrderPtr
        particle_system(thread_, v_0, xi_r, D_r, xi_phi, D_phi, sigma, rho, alpha, pbc_config_, n_, s_);
    BinaryObserverSecondOrderPtr
        binary_observer(thread_, v_0, xi_r, D_r, xi_phi, D_phi, sigma, rho, alpha, pbc_config_, n_, s_, dt, trial);

    Real t = t_0;
    binary_observer.SaveSystemState(system_state_, system_state_size_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(particle_system, system_state_, t, dt);
      // with the linked list technique keep all positions under periodic boundaries
      pbc_config_.ApplyPeriodicBoundaryConditions(system_state_, system_state_size_);
//    thread_->SynchronizeVector(system_state_);
      binary_observer.SaveSystemState(system_state_, system_state_size_, t);
      binary_observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
    } // t
  } // trial

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}