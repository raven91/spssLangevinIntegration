//
// Created by Nikita Kruk on 27.05.20.
//

#ifndef SPSSLANGEVININTEGRATION_PARTICLESYSTEMFIRSTORDERPTR_HPP
#define SPSSLANGEVININTEGRATION_PARTICLESYSTEMFIRSTORDERPTR_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/ThreadSharedMemory.hpp"

class ParticleSystemFirstOrderPtr
{
 public:

  explicit ParticleSystemFirstOrderPtr(ThreadSharedMemory *thread,
                                       Real v_0,
                                       Real xi,
                                       Real sigma,
                                       Real rho,
                                       Real alpha,
                                       Real D_phi,
                                       PeriodicBoundaryConditions &pbc_config,
                                       int number_of_particles, int number_of_state_variables);
  ~ParticleSystemFirstOrderPtr();

  void SetNaturalFrequencies(const std::vector<Real> &natural_frequencies);

  void EvaluateRhs(const Real *const system_state,
                   const std::vector<Real> &k_prev,
                   std::vector<Real> &k_next,
                   Real k_coef,
                   Real dt);

  void EvaluateInteractionsWithGlobalPolarOrderField(const Real *const system_state,
                                                     const std::vector<Real> &k_prev,
                                                     std::vector<Real> &k_next,
                                                     Real k_coef,
                                                     Real dt);
  void EvaluateInteractionsWithAllPairs(const Real *const system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt);
  void EvaluateInteractionsWithLinkedList(const Real *const system_state,
                                          const std::vector<Real> &k_prev,
                                          std::vector<Real> &k_next,
                                          Real k_coef,
                                          Real dt);
  void EvaluateInteractionsWithVerletNeighborList(const Real *const system_state,
                                                  const std::vector<Real> &k_prev,
                                                  std::vector<Real> &k_next,
                                                  Real k_coef,
                                                  Real dt);
  void EvaluateInteractionsWithLinkedListAndVerletNeighborList(const Real *const system_state,
                                                               const std::vector<Real> &k_prev,
                                                               std::vector<Real> &k_next,
                                                               Real k_coef,
                                                               Real dt);

  [[nodiscard]] Real xi() const
  { return xi_; }
  [[nodiscard]] Real rho() const
  { return rho_; }
  [[nodiscard]] Real D_phi() const
  { return D_phi_; }
  [[nodiscard]] Real alpha() const
  { return alpha_; }

  bool uses_verlet_list_;
  void CalculateMaxDisplacement(const std::vector<Real> &total_displacement);
  bool is_higher_order_sde_integration_;
  [[nodiscard]] const std::vector<Real> &GetComplimentaryAlignmentForce() const;

 private:

  ThreadSharedMemory *thread_;
  PeriodicBoundaryConditions &pbc_config_;
  Real v_0_;
  Real xi_;
  Real sigma_;
  Real rho_;
  Real rho_squared_;
  Real alpha_;
  Real cos_alpha_;
  Real sin_alpha_;
  Real D_phi_;
  int n_;
  int s_;

  Real *rk_system_state_;
  std::vector<Real> alignment_force_;
  std::vector<Real> neighborhood_cardinality_;
  std::vector<Real> complementary_alignment_force_; // for higher order SDE integration
  std::vector<Real> natural_frequencies_;

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

};

#endif //SPSSLANGEVININTEGRATION_PARTICLESYSTEMFIRSTORDERPTR_HPP
