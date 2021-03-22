//
// Created by Nikita Kruk on 23.07.20.
//

#include "AveragingObserver.hpp"

#include <algorithm> // std::transform
#include <functional> // std::plus
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip> // std::setprecision

AveragingObserver::AveragingObserver(int number_of_particles, int number_of_state_variables) :
    n_(number_of_particles),
    s_(number_of_state_variables),
    accumulated_system_state_(number_of_particles * number_of_state_variables, 0.0),
    number_of_accumulated_steps_(0)
{

}

AveragingObserver::~AveragingObserver()
{
  accumulated_system_state_.clear();
}

void AveragingObserver::operator()(const std::vector<Real> &system_state, Real t)
{
  std::transform(system_state.begin(),
                 system_state.end(),
                 accumulated_system_state_.begin(),
                 accumulated_system_state_.begin(),
                 std::plus<>{});
  ++number_of_accumulated_steps_;
}

void AveragingObserver::SaveAveragedSystemState(Real v_0,
                                                Real xi,
                                                Real sigma,
                                                Real rho,
                                                Real alpha,
                                                Real D_phi,
                                                int trial)
{
#if defined(__linux__) && (defined(BCS_CLUSTER))
  std::string
      file_folder("/home/nkruk/cpp/spssLangevinIntegration/output/frequency_splitting_bifurcation/coupling_increase/");
#elif defined(__APPLE__)
  std::string
      file_folder("/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/solitary_fraction_vs_n/");
#endif
  std::ostringstream file_name_buffer;
  file_name_buffer << file_folder << "averaged_angular_velocities_v0_" << v_0 << "_xi_" << xi << "_sigma_" << sigma
                   << "_rho_" << rho << "_alpha_" << std::setprecision(10) << alpha << "_Dphi_" << D_phi
                   << "_N_" << n_ << "_" << trial << ".bin";
  std::ofstream file(file_name_buffer.str(), std::ios::out | std::ios::trunc | std::ios::binary);
  assert(file.is_open());

  std::vector<RealOutput> averaged_angular_velocities(accumulated_system_state_.size(), 0.0);
  for (int i = 0; i < n_; ++i)
  {
    averaged_angular_velocities[i] = accumulated_system_state_[s_ * i + 3] / number_of_accumulated_steps_;
  } // i

  file.write((char *) (&averaged_angular_velocities[0]), n_ * sizeof(RealOutput));
  file.close();
}