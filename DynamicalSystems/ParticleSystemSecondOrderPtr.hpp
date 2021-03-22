//
// Created by Nikita Kruk on 14.02.20.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_SELFPROPELLEDPARTICLESYSTEMWITHINERTIAPTR_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_SELFPROPELLEDPARTICLESYSTEMWITHINERTIAPTR_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

class ParticleSystemSecondOrderPtr
{
 public:

  explicit ParticleSystemSecondOrderPtr(ThreadSharedMemory *thread,
                                        Real v_0,
                                        Real xi_r,
                                        Real D_r,
                                        Real xi_phi,
                                        Real D_phi,
                                        Real sigma,
                                        Real rho,
                                        Real alpha,
                                        PeriodicBoundaryConditions &pbc_config,
                                        int number_of_particles, int number_of_state_variables);
  ~ParticleSystemSecondOrderPtr();

  void EvaluateRhs(const Real *const system_state,
                   const std::vector<Real> &k_prev,
                   std::vector<Real> &k_next,
                   Real k_coef,
                   Real dt);
  void EvaluateRhs(const Real *const system_state,
                   std::vector<Real> &derivative,
                   std::vector<std::vector<Real>> &additional_derivative,
                   Real dt);

  void EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt);
  void EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                        std::vector<Real> &derivative,
                                        std::vector<std::vector<Real>> &additional_derivative,
                                        Real dt);
  void EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                          const std::vector<Real> &k_prev,
                                          std::vector<Real> &k_next,
                                          Real k_coef,
                                          Real dt);
  void EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                          std::vector<Real> &derivative,
                                          std::vector<std::vector<Real>> &additional_derivative,
                                          Real dt);
  void EvaluateInteractionsWithVerletNeighborList(const Real *const system_state,
                                                  const std::vector<Real> &k_prev,
                                                  std::vector<Real> &k_next,
                                                  Real k_coef,
                                                  Real dt);
  void EvaluateInteractionsWithVerletNeighborList(const Real *const system_state,
                                                  std::vector<Real> &derivative,
                                                  std::vector<std::vector<Real>> &additional_derivative,
                                                  Real dt);

  Real rho() const
  { return rho_; }
  Real D_r() const
  { return D_r_; }
  Real D_phi() const
  { return D_phi_; }
  Real alpha() const
  { return alpha_; }

  bool last_coefficient_;

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions &pbc_config_;
  Real v_0_;
  Real xi_r_;
  Real D_r_;
  Real xi_phi_;
  Real D_phi_;
  Real sigma_;
  Real rho_;
  Real rho_squared_;
  Real alpha_;
  int n_;
  int s_;

  Real *rk_system_state_;
  std::vector<Real> alignment_force_;
  std::vector<Real> neighborhood_cardinality_;

  // Cell Subdivision Routine
  Real x_size_;
  Real y_size_;
  int num_subcells_x_;
  int num_subcells_y_;
  int *pre_linked_list_;
  std::vector<std::vector<int>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;

  Real r_max_;
  Real r_min_;
  std::vector<std::vector<int>> verlet_list_;
  Real accumulated_displacement_;
  bool should_update_lists_;

  void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y) const;
  void CalculateMaxDisplacement(const std::vector<Real> &derivative, Real dt);

};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_SELFPROPELLEDPARTICLESYSTEMWITHINERTIAPTR_HPP
