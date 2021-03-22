//
// Created by Nikita Kruk on 06.03.18.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP

#include "../Definitions.hpp"
#include "../Parallelization/Thread.hpp"

class PeriodicBoundaryConditions
{
 public:

  explicit PeriodicBoundaryConditions(Thread *thread,
                                      int number_of_particles,
                                      int number_of_state_variables,
                                      Real x_size,
                                      Real y_size);
  ~PeriodicBoundaryConditions();

  void ClassAEffectiveParticleDistance(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) const;

  void ClassBEffectiveParticleDistanceSigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) const;
  void ClassBEffectiveParticleDistanceUnsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) const;

  void ClassCEffectiveParticleDistanceSigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) const;
  void ClassCEffectiveParticleDistanceUnsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) const;

  void ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc) const;
  void ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state);
  void ApplyPeriodicBoundaryConditions(std::vector<RealOutput> &system_state);
  void ApplyPeriodicBoundaryConditions(Real *const system_state, long size);

  Real GetXSize() const;
  Real GetYSize() const;

 private:

  Thread *thread_;
  int n_;
  int s_;
  Real x_size_;
  Real y_size_;
  Real x_rsize_;
  Real y_rsize_;

};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP
