//
// Created by Nikita Kruk on 27.05.20.
//

#include "ParticleSystemFirstOrder.hpp"

#include <cmath>
#include <algorithm> // std::fill
#include <complex>
#include <unordered_map>

ParticleSystemFirstOrder::ParticleSystemFirstOrder(Thread *thread,
                                                   Real v_0,
                                                   Real xi,
                                                   Real sigma,
                                                   Real rho,
                                                   Real alpha,
                                                   Real D_phi,
                                                   PeriodicBoundaryConditions &pbc_config,
                                                   int number_of_particles, int number_of_state_variables) :
    thread_(thread),
    pbc_config_(pbc_config),
    v_0_(v_0),
    xi_(xi),
    sigma_(sigma),
    rho_(rho),
    rho_squared_(rho * rho),
    alpha_(alpha),
    cos_alpha_(std::cos(alpha)),
    sin_alpha_(std::sin(alpha)),
    D_phi_(D_phi),
    alignment_force_(number_of_particles, 0.0),
    neighborhood_cardinality_(number_of_particles, 0.0),
    natural_frequencies_(number_of_particles, 0.0),
    x_size_(1.0),
    y_size_(1.0),
    num_subcells_x_(int(1.0 / rho)),
    num_subcells_y_(int(1.0 / rho)),
    n_(number_of_particles),
    s_(number_of_state_variables),
    uses_verlet_list_(false)
{
  pre_linked_list_ = std::vector<int>(number_of_particles, 0);
  linked_list_ = std::vector<std::vector<int>>(num_subcells_x_ * num_subcells_y_, std::vector<int>());

  //all neighbors
  neighboring_cells_ =
      {
          {-1, -1}, {0, -1}, {1, -1},
          {-1, 0}, {0, 0}, {1, 0},
          {-1, 1}, {0, 1}, {1, 1}
      };
  //half of neighbors
//	neighboring_cells_ =
//	{
//			{0, 0}, {1, 0},
//			{-1, 1}, {0, 1}, {1, 1}
//	};

  // Verlet neighbor list
  r_min_ = rho;
  r_max_ = rho + 0.05;
  should_update_lists_ = true;
  accumulated_displacement_ = 0.0;
  verlet_list_ = std::vector<std::vector<int>>(thread_->GetNumberOfParticlesPerMpichThread(), std::vector<int>());
}

ParticleSystemFirstOrder::~ParticleSystemFirstOrder()
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
}

void ParticleSystemFirstOrder::SetNaturalFrequencies(const std::vector<Real> &natural_frequencies)
{
  std::copy(natural_frequencies.begin(), natural_frequencies.end(), natural_frequencies_.begin());
}

void ParticleSystemFirstOrder::EvaluateRhs(const std::vector<Real> &system_state,
                                           const std::vector<Real> &k_prev,
                                           std::vector<Real> &k_next,
                                           Real k_coef,
                                           Real dt)
{
  if (rho_ >= M_SQRT1_2)
  {
    EvaluateInteractionsWithGlobalPolarOrderField(system_state, k_prev, k_next, k_coef, dt);
  } else if (rho_ <= 0.25)
  {
//    EvaluateInteractionsWithLinkedList(system_state, k_prev, k_next, k_coef, dt);
    EvaluateInteractionsWithLinkedListAndVerletNeighborList(system_state, k_prev, k_next, k_coef, dt);
    uses_verlet_list_ = true;
  } else if (rho_ <= 0.5)
  {
    EvaluateInteractionsWithVerletNeighborList(system_state, k_prev, k_next, k_coef, dt);
    uses_verlet_list_ = true;
  } else
  {
    EvaluateInteractionsWithAllPairs(system_state, k_prev, k_next, k_coef, dt);
  }
}

void ParticleSystemFirstOrder::operator()(const std::vector<Real> &system_state,
                                          std::vector<Real> &derivative,
                                          const Real dt)
{
  static const std::vector<Real> tmp(system_state.size(), 0.0);
  EvaluateRhs(system_state, tmp, derivative, 0.0, dt);
}

void ParticleSystemFirstOrder::EvaluateInteractionsWithGlobalPolarOrderField(const std::vector<Real> &system_state,
                                                                             const std::vector<Real> &k_prev,
                                                                             std::vector<Real> &k_next,
                                                                             Real k_coef,
                                                                             Real dt)
{
  Real /*x_i = 0.0, y_i = 0.0,*/ phi_i = 0.0, omega_i = 0.0;

  // compute global polar order parameter
  std::complex<Real> polar_order_field(0.0, 0.0);
  for (int i = 0; i < n_; ++i)
  {
    phi_i = system_state[s_ * i + 2] + k_coef * dt * k_prev[s_ * i + 2];
    polar_order_field += std::complex<Real>(std::cos(phi_i), std::sin(phi_i));
  } // i
  polar_order_field /= n_;

  // compute interactions using the global polar order parameter
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    int ii = s_ * i;

//    x_i = system_state[ii] + k_coef * dt * k_prev[ii];
//    y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    omega_i = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];

    k_next[ii] = v_0_ * std::cos(phi_i);
    k_next[ii + 1] = v_0_ * std::sin(phi_i);
    k_next[ii + 2] = omega_i;
//    k_next[ii + 3] =
//        -xi_ * omega_i + sigma_ * std::abs(polar_order_field) * std::sin(std::arg(polar_order_field) - phi_i - alpha_);
    k_next[ii + 3] = -xi_ * omega_i + natural_frequencies_[i] + sigma_ * std::abs(polar_order_field)
        * (std::sin(std::arg(polar_order_field) - phi_i) * cos_alpha_
            - std::cos(std::arg(polar_order_field) - phi_i) * sin_alpha_);
  } // i
}

// Calculate the deterministic part of ODE/SDE
void ParticleSystemFirstOrder::EvaluateInteractionsWithAllPairs(const std::vector<Real> &system_state,
                                                                const std::vector<Real> &k_prev,
                                                                std::vector<Real> &k_next,
                                                                Real k_coef,
                                                                Real dt)
{
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0, omega_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0/*, omega_j = 0.0*/;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  static std::vector<Real> rk_system_state(system_state.size(), Real(0.0));
  for (int i : thread_->GetLoopIndices())
  {
    int ii = i * s_;
    // Note system_state is intentionally kept unsynchronized
    rk_system_state[ii] = system_state[ii] + k_coef * dt * k_prev[ii];
    rk_system_state[ii + 1] = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    rk_system_state[ii + 2] = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    rk_system_state[ii + 3] = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];
  } // i

  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    int ii = s_ * i;

    x_i = rk_system_state[ii];
    y_i = rk_system_state[ii + 1];
    phi_i = rk_system_state[ii + 2];
    omega_i = rk_system_state[ii + 3];

    // j ~ i
    for (int j = 0; j < n_; ++j)
    {
      int jj = s_ * j;

      if (j != i)
      {
        x_j = rk_system_state[jj];
        y_j = rk_system_state[jj + 1];
        phi_j = rk_system_state[jj + 2];
//        omega_j = rk_system_state[jj + 3];

        pbc_config_.ClassCEffectiveParticleDistanceUnsigned(x_i, y_i, x_j, y_j, dx, dy);
        dr_squared = dx * dx + dy * dy;

        if (dr_squared <= rho_squared_)
        {
//          alignment_force_[i] += std::sin(phi_j - phi_i - alpha_);
          alignment_force_[i] += std::sin(phi_j - phi_i) * cos_alpha_ - std::cos(phi_j - phi_i) * sin_alpha_;
          ++neighborhood_cardinality_[i];
        }
      }
    } // j

    k_next[ii] = v_0_ * std::cos(phi_i);
    k_next[ii + 1] = v_0_ * std::sin(phi_i);
    k_next[ii + 2] = omega_i;
    if (neighborhood_cardinality_[i] != 0.0)
    {
      k_next[ii + 3] = -xi_ * omega_i + sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
}

// Calculate the deterministic part of ODE/SDE
void ParticleSystemFirstOrder::EvaluateInteractionsWithLinkedList(const std::vector<Real> &system_state,
                                                                  const std::vector<Real> &k_prev,
                                                                  std::vector<Real> &k_next,
                                                                  Real k_coef,
                                                                  Real dt)
{
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0, omega_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0/*, omega_j = 0.0*/;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  static std::vector<Real> rk_system_state(system_state.size(), Real(0.0));
  for (int i : thread_->GetLoopIndices())
  {
    int ii = i * s_;
    // Note system_state is intentionally kept unsynchronized
    rk_system_state[ii] = system_state[ii] + k_coef * dt * k_prev[ii];
    rk_system_state[ii + 1] = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    rk_system_state[ii + 2] = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    rk_system_state[ii + 3] = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];
  } // i

  for (std::vector<int> &linked_list : linked_list_)
  {
    linked_list.clear();
  } // linked_list

  // construct a linked list
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0;
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    ii = s_ * i;
    x_i = rk_system_state[ii];
    y_i = rk_system_state[ii + 1];
//    phi_i = rk_system_state[ii + 2];
//    omega_i = rk_system_state[ii + 3];

//    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//    {
    pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    pre_linked_list_[i] = i_cell;
  } // i
  thread_->SynchronizePrelinkedList(pre_linked_list_);

  for (int i = 0; i < n_; ++i)
  {
    linked_list_[pre_linked_list_[i]].push_back(i);
  } // i

  // loop using the linked list
  int jj = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  for (int i : loop_indices)
  {
    ii = s_ * i;
    x_i = rk_system_state[ii];
    y_i = rk_system_state[ii + 1];
    phi_i = rk_system_state[ii + 2];
    omega_i = rk_system_state[ii + 3];

//    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//    {
    pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    for (auto &neighboring_cell : neighboring_cells_)
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
          x_j = rk_system_state[jj];
          y_j = rk_system_state[jj + 1];
          phi_j = rk_system_state[jj + 2];
//          omega_j = rk_system_state[jj + 3];

//          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
//          {
//            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
//          }

          pbc_config_.ClassCEffectiveParticleDistanceUnsigned(x_i, y_i, x_j, y_j, dx, dy);
          dr_squared = dx * dx + dy * dy;

          if (dr_squared <= rho_squared_)
          {
//            alignment_force_[i] += std::sin(phi_j - phi_i - alpha_);
            alignment_force_[i] += std::sin(phi_j - phi_i) * cos_alpha_ - std::cos(phi_j - phi_i) * sin_alpha_;
            ++neighborhood_cardinality_[i];
          }
        }
      } // j
    } // neighboring_cell

    k_next[ii] = v_0_ * std::cos(phi_i);
    k_next[ii + 1] = v_0_ * std::sin(phi_i);
    k_next[ii + 2] = omega_i;
    if (neighborhood_cardinality_[i] != 0)
    {
      k_next[ii + 3] = -xi_ * omega_i + sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
}

void ParticleSystemFirstOrder::EvaluateInteractionsWithVerletNeighborList(const std::vector<Real> &system_state,
                                                                          const std::vector<Real> &k_prev,
                                                                          std::vector<Real> &k_next,
                                                                          Real k_coef,
                                                                          Real dt)
{
  static std::vector<Real> rk_system_state(system_state.size(), Real(0.0));
  for (int i : thread_->GetLoopIndices())
  {
    int ii = i * s_;
    // Note system_state is intentionally kept unsynchronized
    rk_system_state[ii] = system_state[ii] + k_coef * dt * k_prev[ii];
    rk_system_state[ii + 1] = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    rk_system_state[ii + 2] = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    rk_system_state[ii + 3] = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];
  } // i

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0, omega_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0/*, omega_j = 0.0*/;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  // if it is time to update the Verlet list
  if (should_update_lists_)
  {
    for (auto &verlet_list_entry : verlet_list_)
    {
      verlet_list_entry.clear();
    } // verlet_list_entry

    const std::vector<int> &loop_indices = thread_->GetLoopIndices();
    for (int i : loop_indices)
    {
      int ii = s_ * i;
      x_i = rk_system_state[ii];
      y_i = rk_system_state[ii + 1];

//      if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//      {
//        pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//      }

      // j ~ i
      for (int j = 0; j < n_; ++j)
      {
        int jj = s_ * j;
        if (j != i)
        {
          x_j = rk_system_state[jj];
          y_j = rk_system_state[jj + 1];

//          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
//          {
//            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
//          }

          pbc_config_.ClassCEffectiveParticleDistanceUnsigned(x_i, y_i, x_j, y_j, dx, dy);
          dr_squared = dx * dx + dy * dy;

          if (dr_squared <= r_max_ * r_max_)
          {
            verlet_list_[i].push_back(j);
          }
        }
      } // j
    } // i

    should_update_lists_ = false;
  } // if should_update_lists

  // loop using the Verlet list
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  int ii = 0, jj = 0;
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    ii = s_ * i;
    x_i = rk_system_state[ii];
    y_i = rk_system_state[ii + 1];
    phi_i = rk_system_state[ii + 2];
    omega_i = rk_system_state[ii + 3];

//    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//    {
//      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//    }

    for (int j : verlet_list_[i])
    {
      jj = s_ * j;
      if (j != i)
      {
        x_j = rk_system_state[jj];
        y_j = rk_system_state[jj + 1];
        phi_j = rk_system_state[jj + 2];

//        if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
//        {
//          pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
//        }

        pbc_config_.ClassCEffectiveParticleDistanceUnsigned(x_i, y_i, x_j, y_j, dx, dy);
        dr_squared = dx * dx + dy * dy;

        if (dr_squared <= rho_squared_)
        {
//          alignment_force_[i] += std::sin(phi_j - phi_i - alpha_);
          alignment_force_[i] += std::sin(phi_j - phi_i) * cos_alpha_ - std::cos(phi_j - phi_i) * sin_alpha_;
          ++neighborhood_cardinality_[i];
        }
      }
    } // j

    k_next[ii] = v_0_ * std::cos(phi_i);
    k_next[ii + 1] = v_0_ * std::sin(phi_i);
    k_next[ii + 2] = omega_i;
    if (neighborhood_cardinality_[i] != 0)
    {
      k_next[ii + 3] = -xi_ * omega_i + sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
}

void ParticleSystemFirstOrder::EvaluateInteractionsWithLinkedListAndVerletNeighborList(const std::vector<Real> &system_state,
                                                                                       const std::vector<Real> &k_prev,
                                                                                       std::vector<Real> &k_next,
                                                                                       Real k_coef,
                                                                                       Real dt)
{
  static std::vector<Real> rk_system_state(system_state.size(), Real(0.0));
  for (int i : thread_->GetLoopIndices())
  {
    int ii = i * s_;
    // Note system_state is intentionally kept unsynchronized
    rk_system_state[ii] = system_state[ii] + k_coef * dt * k_prev[ii];
    rk_system_state[ii + 1] = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    rk_system_state[ii + 2] = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    rk_system_state[ii + 3] = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];
  } // i

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0, omega_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0/*, omega_j = 0.0*/;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  // if it is time to update the linked list and the Verlet list
  if (should_update_lists_)
  {
    for (auto &linked_list_entry : linked_list_)
    {
      linked_list_entry.clear();
    } // linked_list_entry

    // construct a linked list
    int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
    int ii = 0;
    const std::vector<int> &loop_indices = thread_->GetLoopIndices();
    for (int i : loop_indices)
    {
      ii = s_ * i;
      x_i = rk_system_state[ii];
      y_i = rk_system_state[ii + 1];

//      if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//      {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//      }

      i_cell_x = int(x_i / x_size_ * num_subcells_x_);
      i_cell_y = int(y_i / y_size_ * num_subcells_y_);
      i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

      pre_linked_list_[i] = i_cell;
    } // i

    // obtain all the linked_list
    for (int i = 0; i < n_; ++i)
    {
      linked_list_[pre_linked_list_[i]].push_back(i);
    } // i

    // construct the Verlet list using the linked list
    for (auto &verlet_list_entry : verlet_list_)
    {
      verlet_list_entry.clear();
    } // verlet_list_entry

    int jj = 0;
    int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
    for (int i : loop_indices)
    {
      ii = s_ * i;
      x_i = rk_system_state[ii];
      y_i = rk_system_state[ii + 1];

//      if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//      {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//      }

      i_cell_x = int(x_i / x_size_ * num_subcells_x_);
      i_cell_y = int(y_i / y_size_ * num_subcells_y_);
      i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

      for (auto &neighboring_cell : neighboring_cells_)
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
            x_j = rk_system_state[jj];
            y_j = rk_system_state[jj + 1];

//            if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
//            {
//              pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
//            }

            pbc_config_.ClassCEffectiveParticleDistanceUnsigned(x_i, y_i, x_j, y_j, dx, dy);
            dr_squared = dx * dx + dy * dy;

            if (dr_squared <= r_max_ * r_max_)
            {
              verlet_list_[i].push_back(j);
            }
          }
        } // j
      } // neighboring_cell
    } // i

    should_update_lists_ = false;
  } // if should_update_lists

  // loop using the Verlet list
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0, jj = 0;
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    ii = s_ * i;
    x_i = rk_system_state[ii];
    y_i = rk_system_state[ii + 1];
    phi_i = rk_system_state[ii + 2];
    omega_i = rk_system_state[ii + 3];

//    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
//    {
//      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
//    }

    for (int j : verlet_list_[i])
    {
      jj = s_ * j;
      if (j != i)
      {
        x_j = rk_system_state[jj];
        y_j = rk_system_state[jj + 1];
        phi_j = rk_system_state[jj + 2];
//        omega_j = rk_system_state[jj + 3];

//        if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
//        {
//          pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
//        }

        pbc_config_.ClassCEffectiveParticleDistanceUnsigned(x_i, y_i, x_j, y_j, dx, dy);
        dr_squared = dx * dx + dy * dy;

        if (dr_squared <= rho_squared_)
        {
//          alignment_force_[i] += std::sin(phi_j - phi_i - alpha_);
          alignment_force_[i] += std::sin(phi_j - phi_i) * cos_alpha_ - std::cos(phi_j - phi_i) * sin_alpha_;
          ++neighborhood_cardinality_[i];
        }
      }
    } // j

    k_next[ii] = v_0_ * std::cos(phi_i);
    k_next[ii + 1] = v_0_ * std::sin(phi_i);
    k_next[ii + 2] = omega_i;
    if (neighborhood_cardinality_[i] != 0)
    {
      k_next[ii + 3] = -xi_ * omega_i + sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
}

void ParticleSystemFirstOrder::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y)
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

void ParticleSystemFirstOrder::CalculateMaxDisplacement(const std::vector<Real> &total_displacement)
{
  Real max2 = 0.0;
  Real dx = 0.0, dy = 0.0, dist2 = 0.0;
  const std::vector<int> &loop_indices = thread_->GetLoopIndices();
  for (int i : loop_indices)
  {
    dx = total_displacement[s_ * i + 0];
    dy = total_displacement[s_ * i + 1];
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
//  int should_update_lists_global = 0;
//  MPI_Allreduce(&should_update_lists_local, &should_update_lists_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//  should_update_lists_ = (should_update_lists_global > 0);
  should_update_lists_ = should_update_lists_local;
}

void ParticleSystemFirstOrder::AddNoise(const std::vector<Real> &system_state,
                                        std::vector<Real> &derivative,
                                        Real t) const
{
  static std::normal_distribution<Real> norm_dist(0.0, 1.0);
  static const Real intensity_phi = std::sqrt(2.0 * D_phi_);
  for (int i : thread_->GetLoopIndices())
  {
    derivative[s_ * i + 3] = intensity_phi * norm_dist(mersenne_twister_generator);
  } // i
}