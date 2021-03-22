//
// Created by Nikita Kruk on 14.02.20.
//

#include "ParticleSystemSecondOrderPtr.hpp"

#include <cmath>
#include <algorithm> // std::fill
#include <string>

ParticleSystemSecondOrderPtr::ParticleSystemSecondOrderPtr(ThreadSharedMemory *thread,
                                                           Real v_0,
                                                           Real xi_r,
                                                           Real D_r,
                                                           Real xi_phi,
                                                           Real D_phi,
                                                           Real sigma,
                                                           Real rho,
                                                           Real alpha,
                                                           PeriodicBoundaryConditions &pbc_config,
                                                           int number_of_particles, int number_of_state_variables) :
    thread_(thread),
    pbc_config_(pbc_config),
    v_0_(v_0),
    xi_r_(xi_r),
    D_r_(D_r),
    xi_phi_(xi_phi),
    D_phi_(D_phi),
    sigma_(sigma),
    rho_(rho),
    rho_squared_(rho * rho),
    alpha_(alpha),
    alignment_force_(number_of_particles, 0.0),
    neighborhood_cardinality_(number_of_particles, 0.0),
    x_size_(1.0),
    y_size_(1.0),
    num_subcells_x_(int(1.0 / rho)),
    num_subcells_y_(int(1.0 / rho)),
    n_(number_of_particles),
    s_(number_of_state_variables)
{
  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(number_of_state_variables * number_of_particles,
                                  rk_system_state_,
                                  std::string("rk_system_state_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread_->AllocateSharedWindow(0, rk_system_state_, std::string("rk_system_state_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("rk_system_state_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &rk_system_state_);
  }

//  pre_linked_list_ = std::vector<int>(kN, 0);
  if (thread_->IsSharedRoot())
  {
    thread_->AllocateSharedWindow(number_of_particles, pre_linked_list_, std::string("pre_linked_list_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
  } else
  {
    thread_->AllocateSharedWindow(0, pre_linked_list_, std::string("pre_linked_list_window"));
    MPI_Barrier(thread_->GetSharedCommunicator());
    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(thread_->GetWindow(std::string("pre_linked_list_window")),
                         thread_->GetSharedRootRank(),
                         &size,
                         &disp_unit,
                         &pre_linked_list_);
  }
  linked_list_ = std::vector<std::vector<int>>(num_subcells_x_ * num_subcells_y_, std::vector<int>());

  // all neighbors
  neighboring_cells_ =
      {
          {-1, -1}, {0, -1}, {1, -1},
          {-1, 0}, {0, 0}, {1, 0},
          {-1, 1}, {0, 1}, {1, 1}
      };
  // half of neighbors
//	neighboring_cells_ =
//	{
//			{0, 0}, {1, 0},
//			{-1, 1}, {0, 1}, {1, 1}
//	};

  //Verlet neighbor list
  last_coefficient_ = true;
  r_min_ = rho;
  r_max_ = 2 * rho;
  should_update_lists_ = true;
  accumulated_displacement_ = 0.0;
  // the verlet list is not required to be synchronized throughout threads
  verlet_list_ = std::vector<std::vector<int>>(number_of_particles, std::vector<int>());
}

ParticleSystemSecondOrderPtr::~ParticleSystemSecondOrderPtr()
{
  alignment_force_.clear();
  neighborhood_cardinality_.clear();
  neighboring_cells_.clear();
  for (std::vector<int> &linked_list : linked_list_)
  {
    linked_list.clear();
  } // linked_list
  linked_list_.clear();
  for (std::vector<int> &verlet_list : verlet_list_)
  {
    verlet_list.clear();
  } // verlet_list
  verlet_list_.clear();

  MPI_Barrier(thread_->GetSharedCommunicator());
  thread_->FreeSharedWindow(std::string("rk_system_state_window"));
  thread_->FreeSharedWindow(std::string("pre_linked_list_window"));
}

void ParticleSystemSecondOrderPtr::EvaluateRhs(const Real *const system_state,
                                               const std::vector<Real> &k_prev,
                                               std::vector<Real> &k_next,
                                               Real k_coef,
                                               Real dt)
{
  if (rho_ <= 0.25)
  {
    EvaluateInteractionsWithLinkedList(system_state, k_prev, k_next, k_coef, dt);
  } else
  {
    EvaluateInteractionsWithAllPairs(system_state, k_prev, k_next, k_coef, dt);
  }
}

void ParticleSystemSecondOrderPtr::EvaluateRhs(const Real *const system_state,
                                               std::vector<Real> &derivative,
                                               std::vector<std::vector<Real> > &additional_derivative,
                                               Real dt)
{
  if (rho_ <= 0.25)
  {
    EvaluateInteractionsWithLinkedList(system_state, derivative, additional_derivative, dt);
  } else
  {
    EvaluateInteractionsWithAllPairs(system_state, derivative, additional_derivative, dt);
  }
//  EvaluateInteractionsWithVerletNeighborList(system_state, derivative, additional_derivative, dt);
}

void ParticleSystemSecondOrderPtr::EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                                                    const std::vector<Real> &k_prev,
                                                                    std::vector<Real> &k_next,
                                                                    Real k_coef,
                                                                    Real dt)
{
  MPI_Win win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    for (int s = 0; s < s_; ++s)
    {
      rk_system_state_[i * s_ + s] = system_state[i * s_ + s] + k_coef * dt * k_prev[i * s_ + s];
    } // s
  } // i
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, v_x_i = 0.0, v_y_i = 0.0, phi_i = 0.0, omega_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, v_x_j = 0.0, v_y_j = 0.0, phi_j = 0.0, omega_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  MPI_Win_lock_all(0, win_rk_sys);
  for (int i : loop_indices)
  {
    int ii = s_ * i;
    x_i = system_state[ii] + k_coef * dt * k_prev[ii];
    y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    v_x_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    v_y_i = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];
    phi_i = system_state[ii + 4] + k_coef * dt * k_prev[ii + 4];
    omega_i = system_state[ii + 5] + k_coef * dt * k_prev[ii + 5];

    // j ~ i
    for (int j = 0; j < n_; ++j)
    {
      int jj = s_ * j;
      if (j != i)
      {
        x_j = system_state[jj] + k_coef * dt * k_prev[jj];
        y_j = system_state[jj + 1] + k_coef * dt * k_prev[jj + 1];
        v_x_j = system_state[jj + 2] + k_coef * dt * k_prev[jj + 2];
        v_y_j = system_state[jj + 3] + k_coef * dt * k_prev[jj + 3];
        phi_j = system_state[jj + 4] + k_coef * dt * k_prev[jj + 4];
        omega_j = system_state[jj + 5] + k_coef * dt * k_prev[jj + 5];

        pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
        dr_squared = dx * dx + dy * dy;

        if (dr_squared <= rho_squared_)
        {
          alignment_force_[i] += std::sin(phi_j - phi_i - alpha_);
          ++neighborhood_cardinality_[i];
        }
      }
    } // j

    k_next[ii] = v_x_i;
    k_next[ii + 1] = v_y_i;
    k_next[ii + 2] = -xi_r_ * v_x_i + v_0_ * std::cos(phi_i);
    k_next[ii + 3] = -xi_r_ * v_y_i + v_0_ * std::sin(phi_i);
    k_next[ii + 4] = omega_i;
    if (neighborhood_cardinality_[i] != 0.0)
    {
      k_next[ii + 5] = -xi_phi_ * omega_i + sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
  MPI_Win_unlock_all(win_rk_sys);
}

void ParticleSystemSecondOrderPtr::EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                                                    std::vector<Real> &derivative,
                                                                    std::vector<std::vector<Real> > &additional_derivative,
                                                                    Real dt)
{

}

void ParticleSystemSecondOrderPtr::EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                                                      const std::vector<Real> &k_prev,
                                                                      std::vector<Real> &k_next,
                                                                      Real k_coef,
                                                                      Real dt)
{
  MPI_Win win_rk_sys = thread_->GetWindow(std::string("rk_system_state_window"));
  MPI_Win_lock_all(0, win_rk_sys);
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    // Note system_state is intentionally kept unsynchronized
    for (int s = 0; s < s_; ++s)
    {
      rk_system_state_[i * s_ + s] = system_state[i * s_ + s] + k_coef * dt * k_prev[i * s_ + s];
    } // s
  } // i
  MPI_Win_sync(win_rk_sys);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(rk_system_state_);

  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, v_x_i = 0.0, v_y_i = 0.0, phi_i = 0.0, omega_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, v_x_j = 0.0, v_y_j = 0.0, phi_j = 0.0, omega_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  for (std::vector<int> &linked_list : linked_list_)
  {
    linked_list.clear();
  } // linked_list

  // construct a linked list
  MPI_Win win_pre_linked_list = thread_->GetWindow(std::string("pre_linked_list_window"));
  MPI_Win_lock_all(0, win_pre_linked_list);
  MPI_Win_lock_all(0, win_rk_sys);
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0;
//	if (first_coefficient_)
//	{
  for (int i : loop_indices)
  {
    ii = s_ * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    v_x_i = rk_system_state_[ii + 2];
    v_y_i = rk_system_state_[ii + 3];
    phi_i = rk_system_state_[ii + 4];
    omega_i = rk_system_state_[ii + 5];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    pre_linked_list_[i] = i_cell;
  } // i
  MPI_Win_sync(win_pre_linked_list);
  MPI_Barrier(thread_->GetSharedCommunicator());
  MPI_Win_unlock_all(win_rk_sys);
  thread_->SynchronizeVectorThroughoutClusters(pre_linked_list_);
  // obtain all the linked_list
  for (int i = 0; i < n_; ++i)
  {
    linked_list_[pre_linked_list_[i]].push_back(i);
  } // i
  MPI_Win_unlock_all(win_pre_linked_list);

  // loop using the linked list
  MPI_Win_lock_all(0, win_rk_sys);
  int jj = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  for (int i : loop_indices)
  {
    ii = s_ * i;
    x_i = rk_system_state_[ii];
    y_i = rk_system_state_[ii + 1];
    v_x_i = rk_system_state_[ii + 2];
    v_y_i = rk_system_state_[ii + 3];
    phi_i = rk_system_state_[ii + 4];
    omega_i = rk_system_state_[ii + 5];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    for (const std::vector<int> &neighboring_cell : neighboring_cells_)
    {
      j_cell_x = i_cell_x + neighboring_cell[0];
      j_cell_y = i_cell_y + neighboring_cell[1];
      AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
      j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

      for (int j : linked_list_[j_cell])
      {
        jj = s_ * j;
        if (j != i)
        {
          x_j = rk_system_state_[jj];
          y_j = rk_system_state_[jj + 1];
          v_x_j = rk_system_state_[jj + 2];
          v_y_j = rk_system_state_[jj + 3];
          phi_j = rk_system_state_[jj + 4];
          omega_j = rk_system_state_[jj + 5];

          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
          {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
          }
          pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
          dr_squared = dx * dx + dy * dy;

          if (dr_squared <= rho_squared_)
          {
            alignment_force_[i] += std::sin(phi_j - phi_i - alpha_);
            ++neighborhood_cardinality_[i];
          }
        }
      } // j
    } // neighboring_cell

    k_next[ii] = v_x_i;
    k_next[ii + 1] = v_y_i;
    k_next[ii + 2] = -xi_r_ * v_x_i + v_0_ * std::cos(phi_i);
    k_next[ii + 3] = -xi_r_ * v_y_i + v_0_ * std::sin(phi_i);
    k_next[ii + 4] = omega_i;
    if (neighborhood_cardinality_[i] != 0)
    {
      k_next[ii + 5] = -xi_phi_ * omega_i + sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
  MPI_Win_unlock_all(win_rk_sys);
}

void ParticleSystemSecondOrderPtr::EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                                                      std::vector<Real> &derivative,
                                                                      std::vector<std::vector<Real> > &additional_derivative,
                                                                      Real dt)
{

}

void ParticleSystemSecondOrderPtr::EvaluateInteractionsWithVerletNeighborList(const Real *const system_state,
                                                                              const std::vector<Real> &k_prev,
                                                                              std::vector<Real> &k_next,
                                                                              Real k_coef,
                                                                              Real dt)
{
  // TODO: to implement the function
}

void ParticleSystemSecondOrderPtr::EvaluateInteractionsWithVerletNeighborList(const Real *const system_state,
                                                                              std::vector<Real> &derivative,
                                                                              std::vector<std::vector<Real> > &additional_derivative,
                                                                              Real dt)
{

}

void ParticleSystemSecondOrderPtr::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y) const
{
  if (cell_x < 0)
  {
    cell_x = num_subcells_x_ - 1;
  } else if (cell_x >= num_subcells_x_)
  {
    cell_x = 0;
  }

  if (cell_y < 0)
  {
    cell_y = num_subcells_y_ - 1;
  } else if (cell_y >= num_subcells_y_)
  {
    cell_y = 0;
  }
}

void ParticleSystemSecondOrderPtr::CalculateMaxDisplacement(const std::vector<Real> &derivative, Real dt)
{
  Real max2 = 0.0;
  Real dx = 0.0, dy = 0.0, dist2 = 0.0;
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    dx = derivative[s_ * i + 0] * dt;
    dy = derivative[s_ * i + 1] * dt;
    dist2 = dx * dx + dy * dy;
    if (max2 < dist2)
    {
      max2 = dist2;
    }
  } // i
  accumulated_displacement_ += std::sqrt(max2);
  int should_update_lists_local = 0;
  if (accumulated_displacement_ > 0.5 * (r_max_ - r_min_))
  {
    accumulated_displacement_ = 0.0;
    should_update_lists_local = 1;
  }
  int should_update_lists_global = 0;
  MPI_Allreduce(&should_update_lists_local, &should_update_lists_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  should_update_lists_ = (should_update_lists_global > 0);
}