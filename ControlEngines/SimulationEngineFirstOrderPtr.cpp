//
// Created by Nikita Kruk on 27.05.20.
//

#include "SimulationEngineFirstOrderPtr.hpp"
#include "../Steppers/RungeKutta4StepperPtr.hpp"
#include "../Steppers/StochasticRungeKuttaStepperPtr.hpp"
#include "../DynamicalSystems/ParticleSystemFirstOrderPtr.hpp"
#include "../Observers/BinaryObserverFirstOrderPtr.hpp"
#include "../Observers/AveragingObserverPtr.hpp"

#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm> // std::copy
#include <iostream>
#include <complex>

SimulationEngineFirstOrderPtr::SimulationEngineFirstOrderPtr(ThreadSharedMemory *thread,
                                                             int number_of_particles,
                                                             int number_of_state_variables) :
    thread_(thread),
    pbc_config_(thread, number_of_particles, number_of_state_variables, 1.0, 1.0),
    system_state_size_(number_of_state_variables * number_of_particles),
    natural_frequencies_(number_of_particles, 0.0),
    n_(number_of_particles),
    s_(number_of_state_variables)
{
  system_state_ = new Real[system_state_size_];
}

SimulationEngineFirstOrderPtr::~SimulationEngineFirstOrderPtr()
{
  delete[] system_state_;
}

void SimulationEngineFirstOrderPtr::InitializeRandomSystemState(Real sigma, Real xi)
{
  if (thread_->IsRoot())
  {
    const Real two_pi = 2.0 * M_PI;

    std::uniform_real_distribution<Real> unif_real_dist_01(0.0, 1.0);
    std::uniform_real_distribution<Real> unif_real_dist_02pi(0.0, two_pi);
    std::uniform_real_distribution<Real> unif_real_dist_frequency(-sigma / xi, sigma / xi);

    for (int i = 0; i < n_; ++i)
    {
      system_state_[s_ * i] = unif_real_dist_01(mersenne_twister_generator);
      system_state_[s_ * i + 1] = unif_real_dist_01(mersenne_twister_generator);
      system_state_[s_ * i + 2] = unif_real_dist_02pi(mersenne_twister_generator);
      system_state_[s_ * i + 3] = unif_real_dist_frequency(mersenne_twister_generator);
    } // i

    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_, system_state_size_);

  if (thread_->IsRoot())
  {
    std::normal_distribution<Real> natural_frequency_dist(0.0, 0.1);
    for (int i = 0; i < n_; ++i)
    {
      natural_frequencies_[i] = natural_frequency_dist(mersenne_twister_generator);
    } // i
  }
  thread_->BroadcastVector(natural_frequencies_);
}

void SimulationEngineFirstOrderPtr::InitializeSystemStateFromPreviousSolution(const std::string &file_name)
{
  if (thread_->IsRoot())
  {
    std::ifstream prev_solution_file(file_name, std::ios::binary | std::ios::in);
    assert(prev_solution_file.is_open());

    prev_solution_file.seekg(0, std::ios::end);
    std::streampos size = prev_solution_file.tellg();
    int t_max = size / ((1l + s_ * n_) * sizeof(RealOutput));
    prev_solution_file.seekg((1l + s_ * n_) * (t_max - 1) * 1 * sizeof(RealOutput), std::ios::beg);

    RealOutput t = 0.0;
    std::vector<RealOutput> output_system_state(s_ * n_, 0.0);

    prev_solution_file.read((char *) &t, sizeof(RealOutput));
    std::cout << "started from " << t << "s" << std::endl;
    prev_solution_file.read((char *) &output_system_state[0], s_ * n_ * sizeof(RealOutput));
    std::copy(output_system_state.begin(), output_system_state.end(), &system_state_[0]);

    prev_solution_file.close();
    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_, system_state_size_);
}

void SimulationEngineFirstOrderPtr::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

  Real v_0 = 1.0;
  Real xi = 0.1;
  Real sigma = 1.0;
  Real rho = 0.8;
  Real alpha = 0.3;
  Real D_phi = 0.01;

  Real t_0 = 0.0;
  Real t_1 = 1000.0;
  Real dt = 0.01;

  int trial = 0;
//  for (int trial = 0; trial < 100; ++trial)
//  for (int ia = 0; ia < 31; ++ia)
  {
//    alpha = 0.05 * ia;
    InitializeRandomSystemState(sigma, xi);
//    InitializeSystemStateFromPreviousSolution("/Volumes/Kruk/spss/spssLangevinIntegration/spiral_waves/v0_0.01_xi_0.1_sigma_1_rho_0.01_alpha_0.5_Dphi_0.01_N_50000_0_0.bin");
//    RungeKutta4StepperPtr stepper(thread_, n_, s_);
    StochasticRungeKuttaStepperPtr stepper(thread_, n_, s_);
    ParticleSystemFirstOrderPtr particle_system(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_);
//    particle_system.SetNaturalFrequencies(natural_frequencies_);
    BinaryObserverFirstOrderPtr
        binary_observer(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_, dt, trial);

    Real t = t_0;
    binary_observer.SaveSystemState(system_state_, system_state_size_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(particle_system, system_state_, t, dt);
      // with the linked list technique keep all positions under periodic boundaries
////      pbc_config_.ApplyPeriodicBoundaryConditions(system_state_, system_state_size_);
////    thread_->SynchronizeVector(system_state_);
      binary_observer.SaveSystemState(system_state_, system_state_size_, t);
      binary_observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
    } // t

    /*// calculate averaged quantities
    if (thread_->IsRoot())
    {
      std::cout << "starting averaging procedure" << std::endl;
    }
    Real *averaging_system_state = new Real[system_state_size_];
    std::copy(&system_state_[0], &system_state_[system_state_size_], &averaging_system_state[0]);
    AveragingObserverPtr averaging_observer(thread_, n_, s_);
    t = 0.0;
    while (t <= 100.0)
    {
      t += dt;
      stepper.DoStep(particle_system, averaging_system_state, t, dt);
//      thread_->SynchronizeVector(averaging_system_state);
      averaging_observer.operator()(averaging_system_state, system_state_size_, t);
    } // t
    averaging_observer.SaveAveragedSystemState(v_0, xi, sigma, rho, alpha, D_phi, trial);
    delete[] averaging_system_state;*/
  } // trial

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}

void SimulationEngineFirstOrderPtr::RunContinuationMethod()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

  Real v_0 = 1.0;
  Real xi = 0.1;
  Real sigma = 1.0;
  Real rho = 0.8;
  Real alpha = 0.17;
  Real D_phi = 0.0;

  Real t_0 = 0.0;
  Real t_1 = 200.0;
  Real dt = 0.01;

  InitializeSystemStateFromPreviousSolution(
      "/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/two_solitary_clusters/v0_1_xi_0.1_sigma_1_rho_0.8_alpha_0.6_Dphi_0_N_1000_0_0.bin");

  int trial = 0;
//  for (int trial = 0; trial < 16; ++trial)
  for (int ia = 1; ia < 6; ++ia)
  {
    alpha = 0.6 - 0.05 * ia;
    RungeKutta4StepperPtr stepper(thread_, n_, s_);
    ParticleSystemFirstOrderPtr particle_system(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_);
    BinaryObserverFirstOrderPtr
        binary_observer(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_, dt, trial);

    Real t = t_0;
//    binary_observer.SaveSystemState(system_state_, system_state_size_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(particle_system, system_state_, t, dt);
      // with the linked list technique keep all positions under periodic boundaries
//      pbc_config_.ApplyPeriodicBoundaryConditions(system_state_, system_state_size_);
//    thread_->SynchronizeVector(system_state_);
//      binary_observer.SaveSystemState(system_state_, system_state_size_, t);
      binary_observer.SaveSummaryStatistics(system_state_, system_state_size_, t);
    } // t
    binary_observer.SaveSystemState(system_state_, system_state_size_, t);
  } // trial

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}