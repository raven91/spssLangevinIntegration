//
// Created by Nikita Kruk on 27.05.20.
//

#include "SimulationEngineFirstOrder.hpp"
#include "../Steppers/RungeKutta4Stepper.hpp"
#include "../DynamicalSystems/ParticleSystemFirstOrder.hpp"
#include "../Observers/BinaryObserverFirstOrder.hpp"
#include "../Observers/AveragingObserver.hpp"

#include <iostream>
#include <cassert>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <boost/numeric/odeint.hpp>

SimulationEngineFirstOrder::SimulationEngineFirstOrder(Thread *thread,
                                                       int number_of_particles,
                                                       int number_of_state_variables) :
    thread_(thread),
    pbc_config_(thread_, number_of_particles, number_of_state_variables, 1.0, 1.0),
    system_state_(number_of_state_variables * number_of_particles, 0.0),
    natural_frequencies_(number_of_particles, 0.0),
    n_(number_of_particles),
    s_(number_of_state_variables)
{

}

SimulationEngineFirstOrder::~SimulationEngineFirstOrder()
{
  system_state_.clear();
}

void SimulationEngineFirstOrder::InitializeRandomSystemState(Real sigma, Real xi)
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
  thread_->BroadcastVector(system_state_);

  if (thread_->IsRoot())
  {
    std::normal_distribution<Real> natural_frequency_dist(0.0, 0.5);
    for (int i = 0; i < n_; ++i)
    {
      natural_frequencies_[i] = natural_frequency_dist(mersenne_twister_generator);
    } // i
  }
  thread_->BroadcastVector(natural_frequencies_);
}

void SimulationEngineFirstOrder::InitializeSystemStateFromPreviousSolution(const std::string &file_name)
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
    std::vector<RealOutput> input_system_state(s_ * n_, 0.0);

    prev_solution_file.read((char *) &t, sizeof(RealOutput));
    std::cout << "started from " << t << "s" << std::endl;
    prev_solution_file.read((char *) &input_system_state[0], s_ * n_ * sizeof(RealOutput));
    std::copy(input_system_state.begin(), input_system_state.end(), &system_state_[0]);

    prev_solution_file.close();
    std::cout << "system state initialization complete" << std::endl;
  }
  thread_->BroadcastVector(system_state_);
}

#if defined(MPI_PARAMETER_SCAN)
#include <mpi.h>
#endif
void SimulationEngineFirstOrder::RunSimulation()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

#if defined(MPI_PARAMETER_SCAN)
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Real v_0 = 1.0;
  Real xi = 0.1;
  Real sigma = 1.0;
  Real rho = 0.8;
  Real alpha = 0.0 + 0.1 * rank;
  Real D_phi = 0.0;
#else
  Real v_0 = 1.0;
  Real xi = 0.1;
  Real sigma = 5.0;
  Real rho = 0.8;
  Real alpha = 0.0;
  Real D_phi = 0.0;
#endif

  Real t_0 = 0.0;
  Real t_1 = 10000.0;
  Real dt = 0.01;
  /*std::ostringstream angular_velocity_file_name_buffer;
  angular_velocity_file_name_buffer << "/home/nkruk/cpp/spssLangevinIntegration/output/angular_velocity"
                                    << "_v0_" << v_0 << "_xi_" << xi << "_sigma_" << sigma << "_rho_" << rho
                                    << "_alpha_" << alpha << "_Dphi_" << D_phi << "_N_" << n_ << ".txt";
  std::ofstream angular_velocity_file(angular_velocity_file_name_buffer.str(), std::ios::out | std::ios::trunc);
  assert(angular_velocity_file.is_open());*/

  int trial = 0;
//  for (int trial = 0; trial < 100; ++trial)
  {
    InitializeRandomSystemState(sigma, xi);
//    InitializeSystemStateFromPreviousSolution(
//        "/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/v0_1_xi_0.1_sigma_1_rho_0.8_alpha_0.8_Dphi_0_N_10000_0_0.bin");
//    RungeKutta4Stepper stepper(thread_, n_, s_);
    ParticleSystemFirstOrder particle_system(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_);
    particle_system.SetNaturalFrequencies(natural_frequencies_);
    BinaryObserverFirstOrder
        binary_observer(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_, dt, trial);

    /*Real t = t_0;
    binary_observer.SaveSystemState(system_state_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(particle_system, system_state_, t, dt);
      thread_->SynchronizeVector(system_state_);
      if (thread_->IsRoot())
      {
        binary_observer.SaveSystemState(system_state_, t);
        binary_observer.SaveSummaryStatistics(system_state_, t);
      }
    } // t*/
    typedef boost::numeric::odeint::runge_kutta_dopri5<std::vector<Real>> dopri5_type;
    boost::numeric::odeint::integrate_const(boost::numeric::odeint::make_dense_output(1.0e-6, 1.0e-6, dopri5_type()), // 1.0e-10 for AIP:Chaos
                                               particle_system,
                                               system_state_,
                                               t_0,
                                               t_1,
                                               dt,
                                               boost::ref(binary_observer));
//    binary_observer.SaveSystemState(system_state_, t_1);

    /*// calculate averaged quantities
    std::cout << "starting averaging procedure" << std::endl;
    std::vector<Real> averaging_system_state(system_state_);
    AveragingObserver averaging_observer(n_, s_);
    boost::numeric::odeint::integrate_const(boost::numeric::odeint::make_dense_output(1.0e-10, 1.0e-10, dopri5_type()),
                                            particle_system,
                                            averaging_system_state,
                                            0.0,
                                            100.0,
                                            dt,
                                            boost::ref(averaging_observer));
    averaging_observer.SaveAveragedSystemState(v_0, xi, sigma, rho, alpha, D_phi, trial);*/

    /*std::cout << "trial " << trial << std::endl;
    for (int i = 0; i < n_; ++i)
    {
      angular_velocity_file << system_state_[s_ * i + 3] << '\t';
    } // i
    angular_velocity_file << std::endl;*/
#if defined(MPI_PARAMETER_SCAN)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  } // trial

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}

#if defined(MPI_PARAMETER_SCAN)
#include <mpi.h>
#endif
void SimulationEngineFirstOrder::RunContinuationMethod()
{
  if (thread_->IsRoot())
  {
    std::cout << "simulation started with " << thread_->GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

#if defined(MPI_PARAMETER_SCAN)
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Real v_0 = 1.0;
  Real xi = 0.1;
  Real sigma = 1.0;
  Real rho = 0.8;
  Real alpha = 0.0; // + 0.1 * rank;
  Real D_phi = 0.0;
#else
  Real v_0 = 1.0;
  Real xi = 0.1;
  Real sigma = 1.0;
  Real rho = 0.8;
  Real alpha = 0.13;
  Real D_phi = 0.0;
#endif

  Real t_0 = 0.0;
  Real t_1 = 1000.0;
  Real dt = 0.01;

  char buffer[50];
  std::sprintf(buffer, "_%g_", alpha);
  InitializeSystemStateFromPreviousSolution(std::string(
      "/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/one_solitary_particle/N100/first_solitary_bifurcation/v0_1_xi_0.1_sigma_1_rho_0.8_alpha")
                                                + buffer + std::string("Dphi_0_N_100_0_0.bin"));

  int trial = 0;
//  for (int trial = 0; trial < 16; ++trial)
  for (int ia = 1; ia < 11; ++ia)
//  for (int is = 1; is < 100; ++is)
  {
    alpha = 0.13 - 0.00001 * ia;
//    sigma = 0.19 - is * 0.0001;
    std::cout << "alpha=" << alpha << ", sigma=" << sigma << std::endl;
    RungeKutta4Stepper stepper(thread_, n_, s_);
    ParticleSystemFirstOrder particle_system(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_);
    BinaryObserverFirstOrder
        binary_observer(thread_, v_0, xi, sigma, rho, alpha, D_phi, pbc_config_, n_, s_, dt, trial);

    /*Real t = t_0;
//    binary_observer.SaveSystemState(system_state_, t);
    while (t <= t_1)
    {
      t += dt;
      stepper.DoStep(particle_system, system_state_, t, dt);
      thread_->SynchronizeVector(system_state_);
      if (thread_->IsRoot())
      {
//        binary_observer.SaveSystemState(system_state_, t);
        binary_observer.SaveSummaryStatistics(system_state_, t);
      }
    } // t*/
    typedef boost::numeric::odeint::runge_kutta_dopri5<std::vector < Real>>
    dopri5_type;
    boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_dense_output(1.0e-10,
                                                                                         1.0e-10,
                                                                                         dopri5_type()),
                                               particle_system,
                                               system_state_,
                                               t_0,
                                               t_1,
                                               dt);
//    binary_observer.SaveSystemState(system_state_, t_1);

    // calculate averaged quantities
    std::cout << "starting averaging procedure" << std::endl;
    std::vector<Real> averaging_system_state(system_state_);
    AveragingObserver averaging_observer(n_, s_);
    boost::numeric::odeint::integrate_const(boost::numeric::odeint::make_dense_output(1.0e-10, 1.0e-10, dopri5_type()),
                                            particle_system,
                                            averaging_system_state,
                                            0.0,
                                            100.0,
                                            dt,
                                            boost::ref(averaging_observer));
    averaging_observer.SaveAveragedSystemState(v_0, xi, sigma, rho, alpha, D_phi, trial);

    boost::numeric::odeint::integrate_const(boost::numeric::odeint::make_dense_output(1.0e-10, 1.0e-10, dopri5_type()),
                                            particle_system,
                                            system_state_,
                                            t_0,
                                            30.0,
                                            dt,
                                            boost::ref(binary_observer));

//#if defined(MPI_PARAMETER_SCAN)
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif
  }

  if (thread_->IsRoot())
  {
    std::cout << "simulation complete" << std::endl;
  }
}