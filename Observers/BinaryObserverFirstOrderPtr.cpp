//
// Created by Nikita Kruk on 27.05.20.
//

#include "BinaryObserverFirstOrderPtr.hpp"

#include <mpi.h>
#include <mpio.h>
#include <sstream>
#include <cassert>
#include <cmath>
#include <complex>
#include <numeric> // std::accumulate
#include <iostream>
#include <iterator> // std::back_inserter, std::copy
#include <boost/math/constants/constants.hpp>

const Real kTwoPi = boost::math::constants::two_pi<Real>();

BinaryObserverFirstOrderPtr::BinaryObserverFirstOrderPtr(ThreadSharedMemory *thread,
                                                         Real v_0,
                                                         Real xi,
                                                         Real sigma,
                                                         Real rho,
                                                         Real alpha,
                                                         Real D_phi,
                                                         PeriodicBoundaryConditions &pbc_config,
                                                         int number_of_particles,
                                                         int number_of_state_variables,
                                                         Real dt,
                                                         int trial) :
    thread_(thread),
    pbc_config_(pbc_config),
    n_(number_of_particles),
    s_(number_of_state_variables),
    output_time_counter_{0, 0},
    output_time_threshold_{100, 100} // mod 1 - save at every dt
{
  integration_step_timer_ = std::chrono::system_clock::now();

#if defined(__linux__) && defined(LICHTENBERG)
  std::string file_folder("/work/scratch/nk59zoce/cpp/spssLangevinIntegration/rk4/");
#elif defined(__linux__) && (defined(BCS_CLUSTER))
  std::string file_folder("/home/nkruk/cpp/spssLangevinIntegration/output/");
#elif defined(__linux__)
  std::string file_folder("/home/nikita/Documents/spssLangevinIntegration/");
#elif defined(__APPLE__)
//  std::string file_folder("/Volumes/Kruk/spss/spssLangevinIntegration/");
  std::string file_folder("/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/");
#endif

  std::ostringstream simulation_file_name_buffer;
  simulation_file_name_buffer << file_folder << "v0_" << v_0 << "_xi_" << xi << "_sigma_" << sigma
                              << "_rho_" << rho << "_alpha_" << alpha << "_Dphi_" << D_phi
                              << "_N_" << number_of_particles;
  simulation_file_name_buffer << "_" << 0;
  simulation_file_name_buffer << "_" << trial << ".bin";
  simulation_file_name_ = simulation_file_name_buffer.str();
  if (thread_->IsRoot())
  {
    std::remove(simulation_file_name_.c_str());
  }
  simulation_file_.open(simulation_file_name_, std::ios::binary | std::ios::out | std::ios::app);
  assert(simulation_file_.is_open());

  std::ostringstream summary_statistics_file_name_buffer;
  summary_statistics_file_name_buffer << file_folder << "v0_" << v_0 << "_xi_" << xi << "_sigma_" << sigma
                                      << "_rho_" << rho << "_alpha_" << alpha << "_Dphi_" << D_phi
                                      << "_N_" << number_of_particles;
  summary_statistics_file_name_buffer << "_" << 0;
  summary_statistics_file_name_buffer << "_" << trial << ".txt";
  summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
  if (thread_->IsRoot())
  {
    std::remove(summary_statistics_file_name_.c_str());
  }
  summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
  assert(summary_statistics_file_.is_open());

  MPI_Barrier(MPI_COMM_WORLD);
}

BinaryObserverFirstOrderPtr::~BinaryObserverFirstOrderPtr()
{
  if (thread_->IsRoot())
  {
    if (simulation_file_.is_open())
    {
      simulation_file_.close();
    }
    if (summary_statistics_file_.is_open())
    {
      summary_statistics_file_.close();
    }
  }
}

void BinaryObserverFirstOrderPtr::SaveSystemState(Real *const system_state, long size, Real t)
{
  thread_->SynchronizeVector(system_state, size);
  if (thread_->IsRoot())
  {
    std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
    std::cout << "value at t = " << t << " integrated in " << elapsed_seconds.count() << "s" << std::endl;
    integration_step_timer_ = std::chrono::system_clock::now();

//	if (t >= 2000.0-100.0)
    if (!(output_time_counter_[0] % output_time_threshold_[0]))
    {
      RealOutput output_t = t;
      static std::vector<RealOutput> output_system_state(size, 0.0);
      std::copy(&system_state[0], &system_state[size], &output_system_state[0]);
      pbc_config_.ApplyPeriodicBoundaryConditions(output_system_state);
      for (int i = 0; i < n_; ++i)
      {
        output_system_state[s_ * i + 2] -= std::floor(output_system_state[s_ * i + 2] / kTwoPi) * kTwoPi;
      } // i

      simulation_file_.write((char *) (&output_t), sizeof(RealOutput));
      simulation_file_.write((char *) (&output_system_state[0]), s_ * n_ * sizeof(RealOutput));
    }
    ++output_time_counter_[0];
  }
}

void BinaryObserverFirstOrderPtr::SaveSummaryStatistics(const Real *const system_state, long size, Real t)
{
  if (thread_->IsRoot())
  {
    if (!(output_time_counter_[1] % output_time_threshold_[1]))
    {
      float t_float = t;
      std::complex<float> complex_order_parameter(0.0f, 0.0f);
      for (int i = 0; i < n_; ++i)
      {
        Real phi_i = system_state[s_ * i + 2];
        complex_order_parameter += std::complex<float>((float) std::cos(phi_i), (float) std::sin(phi_i));
      } // i
      complex_order_parameter /= n_;

      summary_statistics_file_ << t_float << '\t' << std::abs(complex_order_parameter) << '\t'
                               << std::arg(complex_order_parameter) << std::endl;
    }
    ++output_time_counter_[1];
  }
}