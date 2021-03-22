//
// Created by Nikita Kruk on 27.05.20.
//

#ifndef SPSSLANGEVININTEGRATION_PARTICLESYSTEMFIRSTORDER_HPP
#define SPSSLANGEVININTEGRATION_PARTICLESYSTEMFIRSTORDER_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../Parallelization/Thread.hpp"

class ParticleSystemFirstOrder
{
 public:

  explicit ParticleSystemFirstOrder(Thread *thread,
                                    Real v_0,
                                    Real xi,
                                    Real sigma,
                                    Real rho,
                                    Real alpha,
                                    Real D_phi,
                                    PeriodicBoundaryConditions &pbc_config,
                                    int number_of_particles, int number_of_state_variables);
  ~ParticleSystemFirstOrder();

  void SetNaturalFrequencies(const std::vector<Real> &natural_frequencies);

  void EvaluateRhs(const std::vector<Real> &system_state,
                   const std::vector<Real> &k_prev,
                   std::vector<Real> &k_next,
                   Real k_coef,
                   Real dt);
  void operator()(const std::vector<Real> &system_state, std::vector<Real> &derivative, const Real dt);

  void EvaluateInteractionsWithGlobalPolarOrderField(const std::vector<Real> &system_state,
                                                     const std::vector<Real> &k_prev,
                                                     std::vector<Real> &k_next,
                                                     Real k_coef,
                                                     Real dt);
  void EvaluateInteractionsWithAllPairs(const std::vector<Real> &system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt);
  void EvaluateInteractionsWithLinkedList(const std::vector<Real> &system_state,
                                          const std::vector<Real> &k_prev,
                                          std::vector<Real> &k_next,
                                          Real k_coef,
                                          Real dt);
  void EvaluateInteractionsWithVerletNeighborList(const std::vector<Real> &system_state,
                                                  const std::vector<Real> &k_prev,
                                                  std::vector<Real> &k_next,
                                                  Real k_coef,
                                                  Real dt);
  void EvaluateInteractionsWithLinkedListAndVerletNeighborList(const std::vector<Real> &system_state,
                                                               const std::vector<Real> &k_prev,
                                                               std::vector<Real> &k_next,
                                                               Real k_coef,
                                                               Real dt);

  // stochastic part of the Stochastic Euler method
  void AddNoise(const std::vector<Real> &system_state, std::vector<Real> &derivative, Real t) const;

  Real rho() const
  { return rho_; }
  Real D_phi() const
  { return D_phi_; }

  bool uses_verlet_list_;
  void CalculateMaxDisplacement(const std::vector<Real> &total_displacement);

 private:

  Thread *thread_;
  PeriodicBoundaryConditions &pbc_config_;
  Real v_0_;
  Real xi_;
  Real D_phi_;
  Real sigma_;
  Real rho_;
  Real rho_squared_;
  Real alpha_;
  Real cos_alpha_;
  Real sin_alpha_;
  int n_;
  int s_;

  std::vector<Real> alignment_force_;
  std::vector<Real> neighborhood_cardinality_;
  std::vector<Real> natural_frequencies_;

  // Cell Subdivision Routine
  Real x_size_;
  Real y_size_;
  int num_subcells_x_;
  int num_subcells_y_;
  std::vector<int> pre_linked_list_;
  std::vector<std::vector<int>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;

  Real r_max_;
  Real r_min_;
  std::vector<std::vector<int>> verlet_list_;
  Real accumulated_displacement_;
  bool should_update_lists_;

  void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y);

};

#endif //SPSSLANGEVININTEGRATION_PARTICLESYSTEMFIRSTORDER_HPP
