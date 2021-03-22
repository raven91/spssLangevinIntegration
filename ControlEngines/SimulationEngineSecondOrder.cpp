//
// Created by Nikita Kruk on 06.03.18.
//

#include "SimulationEngineSecondOrder.hpp"
#include "../Steppers/RungeKutta4Stepper.hpp"
#include "../DynamicalSystems/ParticleSystemSecondOrder.hpp"
#include "../Observers/BinaryObserverSecondOrder.hpp"

#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm> // std::copy
#include <iostream>

SimulationEngineSecondOrder::SimulationEngineSecondOrder(Thread *thread,
                                                         int number_of_particles,
                                                         int number_of_state_variables) :
    thread_(thread),
    pbc_config_(thread_, number_of_particles, number_of_state_variables, 1.0, 1.0),
    system_state_(number_of_state_variables * number_of_particles, 0.0),
    n_(number_of_particles),
    s_(number_of_state_variables)
{

}

SimulationEngineSecondOrder::~SimulationEngineSecondOrder()
{
  system_state_.clear();
}

void SimulationEngineSecondOrder::InitializeRandomSystemState()
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
  thread_->BroadcastVector(system_state_);
}

void SimulationEngineSecondOrder::InitializeSystemStateFromFile(const std::string &file_name)
{
  if (thread_->IsRoot())
  {
    std::ifstream init_cond_file(file_name, std::ios::binary | std::ios::in);
    assert(init_cond_file.is_open());

    //	std::streampos size = init_cond_file.tellg();
    init_cond_file.seekg(0, std::ios::beg);

//		float t = 0.0f;
    std::vector<float> system_state_float(s_ * n_, 0.0f);

//		init_cond_file.read((char *)&t, sizeof(float));
    init_cond_file.read((char *) &system_state_float[0], s_ * n_ * sizeof(float));
    std::copy(system_state_float.begin(), system_state_float.end(), system_state_.begin());

    init_cond_file.close();
    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_);
}

#if defined(MPI_PARAMETER_SCAN)
#include <mpi.h>
#endif
void SimulationEngineSecondOrder::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

#if defined(MPI_PARAMETER_SCAN)
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Real v_0 = 1.0;
  Real xi_r = 1.0;
  Real D_r = 0.0;
  Real xi_phi = 1.0;
  Real D_phi = 0.0;
  Real sigma = 1.0;
  Real rho = 0.3;
  Real alpha = 1.54;
#else
  Real v_0 = 1.0;
  Real xi_r = 1.0;
  Real D_r = 0.0;
  Real xi_phi = 1.0;
  Real D_phi = 0.0;
  Real sigma = 1.0;
  Real rho = 0.3;
  Real alpha = 1.54;
#endif

  Real t_0 = 0.0;
  Real t_1 = 100.0;
  Real dt = 0.01;
//	Real abs_err = 1.0e-10;
//	Real rel_err = 1.0e-6;

  int trial = 0;
//  for (int trial = 0; trial < 10; ++trial)
  {
    InitializeRandomSystemState();
    RungeKutta4Stepper stepper(thread_, n_, s_);
    ParticleSystemSecondOrder
        particle_system(thread_, v_0, xi_r, D_r, xi_phi, D_phi, sigma, rho, alpha, pbc_config_, n_, s_);
    BinaryObserverSecondOrder
        binary_observer(thread_, v_0, xi_r, D_r, xi_phi, D_phi, sigma, rho, alpha, pbc_config_, n_, s_, dt, trial);

    Real t = t_0;
    binary_observer.SaveSystemState(system_state_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(particle_system, system_state_, t, dt);
      // with the linked list technique keep all positions under periodic boundaries
      pbc_config_.ApplyPeriodicBoundaryConditions(system_state_);
      thread_->SynchronizeVector(system_state_);
      if (thread_->IsRoot())
      {
        binary_observer.SaveSystemState(system_state_, t);
        binary_observer.SaveSummaryStatistics(system_state_, t);
      }
    } // t

#if defined(MPI_PARAMETER_SCAN)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}