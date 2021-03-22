//
// Created by Nikita Kruk on 01.09.20.
//

#include "AveragingObserverPtr.hpp"

#include <algorithm> // std::transform
#include <functional> // std::plus
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip> // std::setprecision

AveragingObserverPtr::AveragingObserverPtr(ThreadSharedMemory *thread,
                                           int number_of_particles,
                                           int number_of_state_variables) :
    thread_(thread),
    n_(number_of_particles),
    s_(number_of_state_variables),
    accumulated_system_state_(number_of_particles * number_of_state_variables, 0.0),
    number_of_accumulated_steps_(0)
{

}

AveragingObserverPtr::~AveragingObserverPtr()
{
  accumulated_system_state_.clear();
}

void AveragingObserverPtr::operator()(Real *const system_state, long size, Real t)
{
  thread_->SynchronizeVector(system_state, size);
  if (thread_->IsRoot())
  {
    std::transform(&system_state[0],
                   &system_state[size],
                   accumulated_system_state_.begin(),
                   accumulated_system_state_.begin(),
                   std::plus<>{});
    ++number_of_accumulated_steps_;
  }
}

void AveragingObserverPtr::SaveAveragedSystemState(Real v_0,
                                                   Real xi,
                                                   Real sigma,
                                                   Real rho,
                                                   Real alpha,
                                                   Real D_phi,
                                                   int trial)
{
  if (thread_->IsRoot())
  {
#if defined(__linux__) && (defined(BCS_CLUSTER))
    std::string
      file_folder("/home/nkruk/cpp/spssLangevinIntegration/output/solitary_fraction_vs_n/n_10000/");
#elif defined(__APPLE__)
    std::string
        file_folder("/Users/nikita/Documents/Projects/spss/spssLangevinIntegration/solitary_fraction_vs_n/n_100/");
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
}