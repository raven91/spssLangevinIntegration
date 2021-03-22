//
// Created by Nikita Kruk on 06.03.18.
//

#include "BinaryObserverSecondOrder.hpp"

#include <sstream>
#include <cassert>
#include <cmath>
#include <complex>
#include <numeric> // std::accumulate
#include <iostream>
#include <iterator> // std::copy

BinaryObserverSecondOrder::BinaryObserverSecondOrder(Thread *thread,
                                                     Real v_0,
                                                     Real xi_r,
                                                     Real D_r,
                                                     Real xi_phi,
                                                     Real D_phi,
                                                     Real sigma,
                                                     Real rho,
                                                     Real alpha,
                                                     PeriodicBoundaryConditions &pbc_config,
                                                     int number_of_particles, int number_of_state_variables,
                                                     Real dt,
                                                     int trial) :
    thread_(thread),
    pbc_config_(pbc_config),
    n_(number_of_particles),
    s_(number_of_state_variables),
    output_time_counter_{0, 0},
    output_time_threshold_{10, 10} // mod 1 - save at every dt
{
  if (thread_->IsRoot())
  {

    integration_step_timer_ = std::chrono::system_clock::now();

#if defined(__linux__) && defined(LICHTENBERG)
    std::string file_folder("/work/scratch/nk59zoce/cpp/spssLangevinIntegration/rk4/");
#elif defined(__linux__) && (defined(BCS_CLUSTER))
    std::string file_folder("/home/nkruk/cpp/spssLangevinIntegration/output/");
#elif defined(__linux__)
    std::string file_folder("/home/nikita/Documents/spssLangevinIntegration/");
#elif defined(__APPLE__)
    std::string file_folder("/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/");
#endif

    std::ostringstream simulation_file_name_buffer;
    simulation_file_name_buffer << file_folder << "v0_" << v_0 << "_xir_" << xi_r
                                << "_Dr_" << D_r << "_xiphi_" << xi_phi << "_Dphi_" << D_phi << "_sigma_"
                                << sigma << "_rho_" << rho << "_alpha_" << alpha << "_N_" << number_of_particles;
//#if defined(MPI_PARAMETER_SCAN)
//    simulation_file_name_buffer << "_" << mpi_thread_rank;
//#else
    simulation_file_name_buffer << "_" << 0;
//#endif
    simulation_file_name_buffer << "_" << trial << ".bin";
    simulation_file_name_ = simulation_file_name_buffer.str();
    std::remove(simulation_file_name_.c_str());

    simulation_file_.open(simulation_file_name_, std::ios::binary | std::ios::out | std::ios::app);
    assert(simulation_file_.is_open());

    std::ostringstream summary_statistics_file_name_buffer;
    summary_statistics_file_name_buffer << file_folder << "summary_statistics_" << "v0_" << v_0 << "_xir_" << xi_r
                                        << "_Dr_" << D_r << "_xiphi_" << xi_phi << "_Dphi_" << D_phi << "_sigma_"
                                        << sigma << "_rho_" << rho << "_alpha_" << alpha << "_N_"
                                        << number_of_particles;
//#if defined(MPI_PARAMETER_SCAN)
//    summary_statistics_file_name_buffer << "_" << mpi_thread_rank;
//#else
    summary_statistics_file_name_buffer << "_" << 0;
//#endif
    summary_statistics_file_name_buffer << "_" << trial << ".txt";
    summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
    std::remove(summary_statistics_file_name_.c_str());

    summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
    assert(summary_statistics_file_.is_open());

  }
}

BinaryObserverSecondOrder::~BinaryObserverSecondOrder()
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

void BinaryObserverSecondOrder::SaveSystemState(std::vector<Real> &system_state, Real t)
{
//	if (t >= 9900.0)
  if (!(output_time_counter_[0] % output_time_threshold_[0]))
  {
    std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
    std::cout << "value at t = " << t << " integrated in " << elapsed_seconds.count() << "s" << std::endl;
    integration_step_timer_ = std::chrono::system_clock::now();

    float t_float = t;
    static std::vector<float> system_state_float(system_state.size(), 0.0);//(system_state.begin(), system_state.end());
    std::copy(system_state.begin(), system_state.end(), system_state_float.begin());

    simulation_file_.write((char *) (&t_float), sizeof(float));
    simulation_file_.write((char *) (&system_state_float[0]), s_ * n_ * sizeof(float));
  }
  ++output_time_counter_[0];
}

void BinaryObserverSecondOrder::SaveSummaryStatistics(std::vector<Real> &system_state, Real t)
{
  if (!(output_time_counter_[1] % output_time_threshold_[1]))
  {
    float t_float = t;
    std::complex<float> complex_order_parameter(0.0f, 0.0f);
    for (int i = 0; i < n_; ++i)
    {
      Real phi_i = system_state[s_ * i + 4];
      complex_order_parameter += std::complex<float>((float) std::cos(phi_i), (float) std::sin(phi_i));
    } // i
    complex_order_parameter /= n_;

    summary_statistics_file_ << t_float << '\t' << std::abs(complex_order_parameter) << '\t'
                             << std::arg(complex_order_parameter) << std::endl;
  }
  ++output_time_counter_[1];
}